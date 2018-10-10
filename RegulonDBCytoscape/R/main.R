#library(dplyr)
#library(RCy3)
#library(plyr)
#library(stringr)



#x es la ruta de la carpeta con los archivos
#x <- c("AgaR")

#x <- c("/Users/daniel/OneDrive\ -\ UNIVERSIDAD\ NACIONAL\ AUTÓNOMA\ DE\ MÉXICO/owncloud/PROYECTO_CYTOSCAPE/Proyecto-Cytoscape/SRsGUs/LacI")

#'@title CytoscapeRegulon
#'@details Creates a cytoscape network from a tab delimited GENSOR unit files.
#'
#'
#' @param x Is the folder path of the GENSOR unit
#' @examples
#' cytoscapeRegulon("/Users/daniel/Proyecto-Cytoscape/SRsGUs/LacI")

#'
#'@author Antonio Daniel Martinez Gutierrez

#' @export
#cytoscapeRegulon <- function(x,html){
cytoscapeRegulon <- function(x){
  #x es el folder con los archivos
  setwd(x)

  ###Obtener el nombre del factor de transcripcion

  nombre_factor_transcripcion <- basename(x)

  reactants_products <- read.delim("reactants_products.txt", header = FALSE,
                                   col.names = c("id","name","type"),
                                   stringsAsFactors=FALSE)


  reactions <- read.delim("reactions.txt", header = FALSE,
                          col.names = c("id","reaction_type"
                                        ,"direction",
                                        "ecocyc_name_id"),
                          stringsAsFactors=FALSE)

  df_nodes <- read.delim("objects.txt", header = FALSE,
                         col.names = c("id","reactant_type"
                                       ,"name"), stringsAsFactors=FALSE)




  complexes <- read.delim("complexes.txt", header = FALSE,
                          col.names = c("id","name","type"),
                          stringsAsFactors=FALSE)

  modifications <- read.delim("modification.txt", header = FALSE,
                              col.names = c("reaction","modification_type"
                                            ,"Object","Reference"), stringsAsFactors=FALSE)



  #####

  ###Identificar las reacciones que necesiten modificación, es decir agregar Nodo "invisible"
  ## y agregar los nodos invisibles a df_nodes
  reac4extraNod<- reactions[reactions$reaction_type!="TRANSLATION","id"]
  num_ex_Nod<-length(reac4extraNod)
  id_Nodaux<-paste0("Aux_",seq_len(num_ex_Nod))
  df_nodes<-rbind(df_nodes,data.frame(id=id_Nodaux,reactant_type="AUX",name=""))
  rm (num_ex_Nod)

  ###Identificar las reacciones que necesiten modificacion y los Objetos que actuan sobre ellas
  ## Establece cuales son reactantes y productos

  names(id_Nodaux)<-reac4extraNod
  mod_nod_prod<-id_Nodaux[modifications$reaction]
  mod_id<-paste0("mo_",seq_len(length(mod_nod_prod)))
  mod_product<-data.frame(id=mod_id,name=mod_nod_prod,type="product")

  mod_nod_reac<-modifications$Object
  mod_id<-paste0("mo_",seq_len(length(mod_nod_reac)))
  mod_reac<-data.frame(id=mod_id,name=mod_nod_reac,type="reactant")

  rm(mod_id,mod_nod_prod,mod_nod_reac)

  ###Los nodos invisibles generan modificaciones en las reacciones que se reconocen a continuacion
  ## Se agregan reacciones temporales
  done<-reactions[(which(reactions$id %in% reac4extraNod,arr.ind = T))*-1,]

  reac_Pr<-reactions[(which(reactions$id %in% reac4extraNod,arr.ind = T)),]
  levels(reac_Pr$id)<-c(levels(reac_Pr$id),paste0(reac4extraNod,"_Pr"))
  reac_Pr$id<-paste0(reac4extraNod,"_Pr")

  reac_Re<-reactions[(which(reactions$id %in% reac4extraNod,arr.ind = T)),]
  levels(reac_Re$id)<-c(levels(reac_Re$id),paste0(reac4extraNod,"_Re"))
  reac_Re$id<-paste0(reac4extraNod,"_Re")

  reactions<-rbind.data.frame(done,reac_Pr,reac_Re)
  rm(done,reac_Pr,reac_Re)

  ###Separar matriz reactant_product por producto y reaccion para
  ## hacer una tercer data.frame formateada para RCY3
  ## Ademas los nodos invisibles cambian los reactantes y productos
  ## se realizan las modificaciones

  df_reactant <- subset(reactants_products, type == "reactant")

  temp_react_id<-paste0(df_reactant[which(df_reactant$id %in% reac4extraNod, arr.ind = T),1],"_Re")
  df_reactant[which(df_reactant$id %in% reac4extraNod, arr.ind = T),]$id <- temp_react_id
  df_reactant <- rbind.data.frame(df_reactant,data.frame(id=paste0(reac4extraNod,"_Pr"),name=id_Nodaux,type="reactant" ),mod_reac)

  df_product <- subset(reactants_products, type == "product")

  temp_product_id<-paste0(df_product[which(df_product$id %in% reac4extraNod, arr.ind = T),1],"_Pr")
  df_product[which(df_product$id %in% reac4extraNod, arr.ind = T),]$id <- temp_product_id
  df_product <- rbind(df_product,data.frame(id=paste0(reac4extraNod,"_Re"),name=id_Nodaux,type="product" ),mod_product)

  df_reactant_product <- merge(df_reactant, df_product, by="id")
  df_reactant_product <- df_reactant_product[,c(1,2,4)]
  colnames(df_reactant_product) <- c("id","reactant_name",
                                     "product_name")

  rm(df_reactant,temp_react_id,df_product,temp_product_id)
  ###Añadir datos de reactions en dataframe con edges final

  df_reactant_product$id<-as.character(df_reactant_product$id)
  reactions$id<-as.character(reactions$id)
  df_edges <- left_join(df_reactant_product, reactions, by = "id")

  colnames(df_edges) <- c("reaction_id","source", "target","interaction",
                          "direction","ecocyc_name")
  ##
  df_edges$reaction_id_ <- df_edges$reaction_id
  df_edges <- df_edges[,-1]
  ##
  #Eliminar dataframes que ya no se usan

  rm(  df_reactant_product, reactants_products, reactions,mod_reac,mod_product)


  ####Añadir que proteinas son TF
  tf_name <- x

  ##extraer rows que son tf
  df_nodes_tf <- subset(df_nodes, id == x)
  df_nodes_tf$reactant_type <- rep("PROTEIN_TF",length(rownames(df_nodes_tf)))

  ##extraer rows que no son tf

  df_nodes_sin_tf <- filter(df_nodes, id != x)

  ###unir nodos ya con tf

  df_nodes <- rbind(df_nodes_tf,df_nodes_sin_tf)

  ###eliminar dataframes que ya no se ocupan

  rm(df_nodes_sin_tf,df_nodes_tf)

  ####





  ############Complejos

  ####Verificar si el archivo de complejos esta vacio
  ####FUNCIONA BIEN, SE COMENTARÁ PARA EVITAR DUPLICACION DE
  ####Miembros del complejo

  #if (nrow(complexes) != 0){
    ############
    #Paso1: obtener los complejos
    #complejos <- unique(complexes$id)
    #Paso2: extraer cada id de complejo
    #for (i in complejos){
      #complejo_df <- complexes[complexes$id==i,]
      #nombre_complejo <- i
      #crear data frame para cada integrante del complejo
      #for (j in seq_along(1:length(rownames(complejo_df)))){
        #integrante_complejo <- complejo_df[j,]
        #nombre_integrante_complejo <- as.character(integrante_complejo$type)
        #df_integrante_complejo <- data.frame("id" = paste0(
          #nombre_integrante_complejo," "),
          #"reactant_type" = "COMPLEX_NODE",
          #"name" = i)
        #HACER UN DF_NODES CON COMPLEJOS
        #if (!exists("df_nodes_complex")){
          #df_nodes_complex <- df_integrante_complejo
        #}


        #Si el dataset existe, se une
        #if (exists("df_nodes_complex")){
          #temp_dataset <- df_integrante_complejo
          #dataset<-rbind(dataset, temp_dataset)
          #df_nodes_complex <- rbind(df_nodes_complex,temp_dataset)
          #rm(temp_dataset)
        #}

      #}

    #}

    #El primer row se repite, se eliminara
    #df_nodes_complex <- df_nodes_complex[2:length(rownames(df_nodes_complex)),]
    #df_nodes_complex$reactant_type <- as.character(df_nodes_complex$reactant_type)


    #Unir el df_nodes_complex con el df_nodes
    #df_nodes <- rbind(df_nodes, df_nodes_complex)
    ######Modificar el df de los complejos con los nuevos nombres con espacio
    #nombres_complejos <- as.character(complexes$type)
    #nombres_complejos <- paste0(nombres_complejos, ' ')
    #complexes$type <- nombres_complejos
    ###Cambiar los colnames de complejes para posteriormente hacerlos edges
    #colnames(complexes) <- c("reaction_id_","target","source")
    #complexes$interaction <- rep("COMPLEX_INTERACTION",
                                 #length(rownames(complexes)))
    ##unir complexes con df_edges
    #df_edges <- rbind.fill(df_edges,complexes)
  #} else{

    #print("INFO: El archivo de complejos esta vacio")

  #}

  ####VERSION 2.0 REVISAR SI LOS COMPLEJOS INTERACCIONAN CON EL TF
  ###########
  if (nrow(complexes) != 0){

    index <- grep(nombre_factor_transcripcion, df_nodes$id)

    reactantes_tf <- df_nodes$reactant_type[index]
    reactantes_tf <- paste(reactantes_tf, "TF", sep="-")

    df_nodes[index,2] <- reactantes_tf
  } else{

  }
  ##########


  #############END COMPLEJOS


  ################Putas Arrows!!!
  ###############

  #Seleccionar primero las que vayan de un reactivo a un aux
  arrow_reactivo_aux <- df_edges[-grep("Aux",df_edges$source),]
  arrow_reactivo_aux$direccion <- rep("reactivo_aux",
                                      length(rownames(arrow_reactivo_aux)))
  #Seleccionar las que vayan de un aux a un reactivo
  arrow_aux_reactivo <- df_edges[grep("Aux",df_edges$source),]
  arrow_aux_reactivo$direccion <- rep("aux_reactivo",
                                      length(rownames(arrow_aux_reactivo)))

  ###Unir las dos df

  arrows_edge <- rbind(arrow_reactivo_aux,
                       arrow_aux_reactivo)
  ###Unir las columnas de interaction con direccion


  arrows_edge$reaction_id_direccion <- gsub('[[:digit:]]+', '', arrows_edge$reaction_id_)

  arrows_edge$arrow <- paste(arrows_edge$interaction,
                             arrows_edge$direccion,sep="-")
  ###Unir las columnas de arrows con direction
  arrows_edge$arrow <- paste(arrows_edge$arrow,
                             arrows_edge$direction,sep="-")
  ##Unir las columnas de arrows con la reaction id direccion
  arrows_edge$arrow <- paste(arrows_edge$arrow,
                             arrows_edge$reaction_id_direccion,sep="-")

  df_edges <- arrows_edge

  rm(arrow_aux_reactivo,arrow_reactivo_aux,arrows_edge)

  #####PRUEBA
  colnames(df_nodes)[colnames(df_nodes)=="name"] <- "description"
  #########


  ###########REMOVER EDGES REPETIDAS

  df_edges$union <- paste(df_edges$source,
                          df_edges$target,sep="-")

  df_edges <- df_edges[!duplicated(df_edges$union),]
  ###REM

  #Llamar RCy3
  createNetworkFromDataFrames(df_nodes,df_edges, title=nombre_factor_transcripcion,
                              collection=nombre_factor_transcripcion)



  ################3Agrupar moleculas de los complejos
  ###Funciona bien, se comentara para facilitar que el grupo de los
  #complejos tenga el color correcto, se puede arreglar.
  #if (nrow(complexes) != 0){
  # nodedata <- getTableColumns("node")
  #nodes_complexes <- nodedata[grep("COMPLEX_NODE", nodedata$reactant_type), ]
  #####
  # Seleccionar los complejos, uno por uno
  #for (g in unique(as.character(complexes$reaction_id_))){
  # deltacatnodes <- df_nodes[grep(g, df_nodes$description), "id"]
  #SUIDS <- nodedata[which(nodedata$`shared name` %in% deltacatnodes),]
  #SUIDS <- as.character(SUIDS$SUID)
  ####
  #nodos_suids <- selectNodes(SUIDS, preserve=FALSE)
  ###Seleccionar el primer nodo vecino
  #primer_vecino <- selectFirstNeighbors()
  #primer_vecino <- setdiff(primer_vecino,nodos_suids)
  #Funcion en r base pero con warnings
  #primer_vecino <- nodedata[which(nodedata$SUID %in% primer_vecino),]
  #primer_vecino <- as.character(primer_vecino$`shared name`)

  ###AQUI EL ERROR DE HTML, AL PARECER EN COLLAPSEGROUP
  ##Crear grupo del complejo
  #    if(missing(html)){
  #createGroup(primer_vecino)###este si corre, no debe estar comentado
  #      collapseGroup(primer_vecino)
  #    }else{
  #      if(html){
  #        createGroup(primer_vecino)
  #      }else{
  #        createGroup(primer_vecino)
  #        collapseGroup(primer_vecino)
  #      }
  #    }
  # }


  #} else{

  #print("INFO: El archivo de complejos esta vacio")
  #}

  ###############


  if (nrow(complexes) != 0){
    column <- 'reactant_type'
    values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
                 'RNA',"GENE","COMPLEX","PROTEIN_TF","COMPLEX_NODE", "AUX",
                 "PROTEIN-TF","COMPLEX-TF","COMPLEX_NODE-TF")
  } else {
    column <- 'reactant_type'
    values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
                 'RNA',"GENE","COMPLEX","PROTEIN_TF","COMPLEX_NODE", "AUX")

  }

  ###Tamaño de los labels




  ########Agregar un default para que los complejos agrupados tengan el estilo
  ####FORMA
  setNodeShapeDefault("OCTAGON")
  setNodeColorDefault("#b6bd7b")
  setNodeBorderColorDefault("#000000")
  setNodeBorderWidthDefault(3)

  #############


  #####Tamaño de los nodos
  lockNodeDimensions(TRUE)


  if (nrow(complexes) != 0){
    sizes <- c(50,50,50,50,50,50,50,10,50,50,50)
    setNodeSizeMapping(column, values, sizes,
                       mapping.type = "d")

  } else{
    sizes <- c(rep(50,7),10)
    setNodeSizeMapping(column, values, sizes,
                       mapping.type = "d")

  }


  #######Forma de los nodos


  if (nrow(complexes) != 0){
    shapes <- c ('ELLIPSE', 'ROUND_RECTANGLE',
                 'PARALLELOGRAM',"RECTANGLE","OCTAGON",
                 "ROUND_RECTANGLE","ELLIPSE","TRIANGLE","ROUND_RECTANGLE",
                 "OCTAGON","ELLIPSE")
    setNodeShapeMapping (column, values, shapes)
  } else{
    shapes <- c ('ELLIPSE', 'ROUND_RECTANGLE',
                 'PARALLELOGRAM',"RECTANGLE","OCTAGON",
                 "ROUND_RECTANGLE","ELLIPSE","TRIANGLE")
    setNodeShapeMapping (column, values, shapes)

  }



  #####Color de los nodos

  if (nrow(complexes) != 0){
    colors <- c("#b6bd7b","#b6bd7b","#ffbc00","#ffbc00","#b6bd7b","#4881a6",
                "#b6bd7b","#ffffff","#4881a6","#4881a6","#4881a6")
    setNodeColorMapping (column, values, colors,
                         mapping.type = "d")

  } else{
    colors <- c("#b6bd7b","#b6bd7b","#ffbc00","#ffbc00","#b6bd7b","#4881a6",
                "#b6bd7b","#ffffff")
    setNodeColorMapping (column, values, colors,
                         mapping.type = "d")

  }


  #####Border de los nodos

  colors <- rep("#000000",length(values))
  setNodeBorderColorMapping (column, values, colors,
                             mapping.type = "d")

  #####Width Border de los nodos

  widths <- rep(3,length(values))
  setNodeBorderWidthMapping(column, values, widths,
                            mapping.type = "d")




  ###Arrow shapes

  #getArrowShapes()
  #?setEdgeLineStyleMapping

  column_line <- "arrow"
  values_line <- c("NA-reactivo_aux-NA-mo_","STATE_TRANSITION-reactivo_aux--re_Re",
                   "TRANSCRIPTION-reactivo_aux--re_RE","STATE_TRANSITION-reactivo_aux-RVB-re_RE",
                   "TRANSLATION-reactivo_aux--re","STATE_TRANSITION-reactivo_aux-L2R-re_RE",
                   "TRANSPORT-reactivo_aux-L2R-re_Re","COMPLEX_INTERACTION-reactivo_aux-NA",
                   "STATE_TRANSITION-aux_reactivo--re_Pr","TRANSCRIPTION-aux_reactivo--re_Pr",
                   "STATE_TRANSITION-aux_reactivo-RVB-re_Pr","STATE_TRANSITION-aux_reactivo-L2R-re_Pr",
                   "TRANSPORT-aux_reactivo-L2R_re_Pr")

  values_line_shape <- c("SOLID","SOLID","SOLID","SOLID","MARQUEE_DASH",
                         "SOLID","SOLID","DASH_DOT","SOLID","SOLID","SOLID")

  setEdgeLineStyleMapping(column_line,values_line, values_line_shape)

  ####PUTAS ARROWS!!!!
  ###Por default
  values_arrow_shape <- c("NONE","NONE","NONE","NONE","OPEN_DELTA","NONE",
                          "NONE","NONE","OPEN_DELTA","OPEN_DELTA","CROSS_DELTA","OPEN_DELTA",
                          "OPEN_DELTA")
  setEdgeTargetArrowMapping(column_line,values_line, values_arrow_shape)

  ###Elegir aquellas edges re_RE que necesiten tener arrow open half circle
  #getArrowShapes
  edgedata <- getTableColumns("edge")
  edges_open_half_circle <- (edgedata[grep("re_Re", df_edges$reaction_id_direccion),"name"])


  setEdgeTargetArrowShapeBypass(edges_open_half_circle,"OPEN_HALF_CIRCLE")








  ####Agregar tabla de modificaciones a cytoscape para hacer los colores de las lineas
  ####Duplicar los valores de las reacciones de modifications y
  #### agregar _Re, y _Pr
  modifications_RE <- modifications
  modifications_RE$reaction <- paste(modifications_RE$reaction, "_Re", sep="")

  modifications_Pr <- modifications
  modifications_Pr$reaction <- paste(modifications_Pr$reaction, "_Pr", sep="")

  modifications <- rbind(modifications_RE, modifications_Pr)

  loadTableData(modifications, data.key.column = "reaction",table = "edge",
                table.key.column = "reaction_id_")

  rm(modifications_RE, modifications_Pr)

  #####Colores de las arrows
  setEdgeTargetArrowColorDefault("#000000")

  column_line_color_arrow <- "modification_type"

  values_color_arrow <- c("PHYSICAL_STIMULATION","INHIBITION","CATALYSIS")

  colores_arrow <-c("#48c4dc","#d80c0c","#848484")

  setEdgeTargetArrowColorMapping(column_line_color_arrow,values_color_arrow,colores_arrow,
                                 mapping.type = "d")

  ######COLORES DE LAS LINEAS

  setEdgeColorDefault("#848484")
  column_line_color_edge <- "modification_type"
  values_line_color_edge <- c("PHYSICAL_STIMULATION","INHIBITION")

  line_color_edge <- c("#48c4dc","#d80c0c")


  setEdgeColorMapping(column_line_color_edge,values_line_color_edge,line_color_edge,
                      mapping.type = "d")


  ###REMOVER TEXTO DE AUX

  nodes_aux <- grep("Aux*",df_nodes$id)
  nodes_aux <- df_nodes[nodes_aux,]
  nodes_aux <- as.character(nodes_aux$id)

  new_nodes_aux_names <- rep(" ", length(nodes_aux))
  setNodeLabelBypass(nodes_aux,new_nodes_aux_names)

  ###TAMAÑOS DE LOS LABELS

  setNodeFontSizeDefault("22")

  nodes_simple_molecule <- df_nodes[grep("SIMPLE_MOLECULE", df_nodes$reactant_type),1]
  setNodeFontSizeBypass(nodes_simple_molecule,16)

  if (nrow(complexes) != 0){
    nodes_complex_label <- df_nodes[grep("COMPLEX*", df_nodes$reactant_type),1]
    setNodeFontSizeBypass(nodes_complex_label,16)

  } else{


  }






  # if(!(missing(html))){
  #  if(html){
  #suppressWarnings(RegulonDBCytoscape::exportToHTML())
  # }
  #}
}
