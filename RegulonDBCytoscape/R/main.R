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
                 "PROTEIN-TF","COMPLEX-TF","COMPLEX_NODE-TF","PROTEIN_TF-TF")
  } else {
    column <- 'reactant_type'
    values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
                 'RNA',"GENE","COMPLEX","PROTEIN_TF","COMPLEX_NODE", "PROTEIN_TF-TF","AUX")

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
  ###Version antigua
  #lockNodeDimensions(TRUE)
  #if (nrow(complexes) != 0){
   # sizes <- c(50,50,50,50,50,50,50,10,50,50,50,50)
  #  setNodeSizeMapping(column, values, sizes,
   #                    mapping.type = "d")

  #} else{
  #  sizes <- c(rep(50,8),10)
  #  setNodeSizeMapping(column, values, sizes,
  #                     mapping.type = "d")

  #}

  ###Version nueva
  lockNodeDimensions(FALSE)
  ###### Proteinas y complejos
  ###### height 35 y width 75
  node_table <- getTableColumns(table = "node")

  if (nrow(complexes) != 0){
    protein_node <- grep("PROTEIN", node_table$reactant_type)
    complex_node <- grep("COMPLEX", node_table$reactant_type)

    index_protein_complex <- sort(c(protein_node, complex_node))
    nodes_protein_complexes <- (node_table[index_protein_complex,"name"])

    setNodeHeightBypass(nodes_protein_complexes, 35)
    setNodeWidthBypass(nodes_protein_complexes, 75)

  } else {

    protein_node <- grep("PROTEIN", node_table$reactant_type)

    nodes_protein <- (node_table[protein_node,"name"])

    setNodeHeightBypass(nodes_protein, 35)
    setNodeWidthBypass(nodes_protein, 75)
  }




  #####Genes y RNAS
  ####35 de altura y 85 de ancho.
  genes_node <- grep("GENE", node_table$reactant_type)

  nodes_genes <- (node_table[genes_node,"name"])

  setNodeHeightBypass(nodes_genes, 35)
  setNodeWidthBypass(nodes_genes, 85)

  #######RNAS
  #####35 de altura y 150 ancho
  rna_node <- grep("RNA", node_table$reactant_type)

  nodes_rna <- (node_table[rna_node,"name"])

  setNodeHeightBypass(nodes_rna, 35)
  setNodeWidthBypass(nodes_rna, 150)

  ###Metabolitos
  ### 35 de altura y 105 de ancho
  metabolites_node <- grep("SIMPLE_MOLECULE", node_table$reactant_type)
  nodes_metabolites <- (node_table[metabolites_node,"name"])

  setNodeHeightBypass(nodes_metabolites, 35)
  setNodeWidthBypass(nodes_metabolites, 105)

  ###Auxiliares
  aux_node <- grep("AUX", node_table$reactant_type)
  nodes_aux <- (node_table[aux_node,"name"])

  setNodeHeightBypass(nodes_aux, 10)
  setNodeWidthBypass(nodes_aux, 10)




  #######Forma de los nodos


  if (nrow(complexes) != 0){
    shapes <- c ('ELLIPSE', 'ROUND_RECTANGLE',
                 'PARALLELOGRAM',"RECTANGLE","OCTAGON",
                 "ROUND_RECTANGLE","ELLIPSE","TRIANGLE","ROUND_RECTANGLE",
                 "OCTAGON","ELLIPSE","ROUND_RECTANGLE")
    setNodeShapeMapping (column, values, shapes)
  } else{
    shapes <- c ('ELLIPSE', 'ROUND_RECTANGLE',
                 'PARALLELOGRAM',"RECTANGLE","OCTAGON",
                 "ROUND_RECTANGLE","ELLIPSE","ROUND_RECTANGLE","TRIANGLE")
    setNodeShapeMapping (column, values, shapes)

  }



  #####Color de los nodos

  if (nrow(complexes) != 0){
    colors <- c("#b6bd7b","#b6bd7b","#ffbc00","#ffbc00","#b6bd7b","#4881a6",
                "#b6bd7b","#ffffff","#4881a6","#4881a6","#4881a6","#4881a6")
    setNodeColorMapping (column, values, colors,
                         mapping.type = "d")

  } else{
    colors <- c("#b6bd7b","#b6bd7b","#ffbc00","#ffbc00","#b6bd7b","#4881a6",
                "#b6bd7b","#4881a6","#ffffff")
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


  setEdgeSourceArrowShapeBypass(edges_open_half_circle,"OPEN_HALF_CIRCLE")




  ####Agregar tabla de modificaciones a cytoscape para hacer los colores de las lineas
  ####Duplicar los valores de las reacciones de modifications y
  #### agregar _Re, y _Pr
  modifications_RE <- modifications
  modifications_RE$reaction <- paste(modifications_RE$reaction, "_Re", sep="")

  modifications_Pr <- modifications
  modifications_Pr$reaction <- paste(modifications_Pr$reaction, "_Pr", sep="")

  modifications <- rbind(modifications_RE, modifications_Pr)

  ####Obtener solo las reacciones unicas


  modifications_unique <- modifications[!duplicated(modifications$reaction), ]

  loadTableData(modifications_unique, data.key.column = "reaction",table = "edge",
                table.key.column = "reaction_id_")

  rm(modifications_RE, modifications_Pr)

  #####Colores de las arrows
  ####MAL, HACERLO EN LA EXTRACCION DE EDGES
  setEdgeTargetArrowColorDefault("#000000")

  column_line_color_arrow <- "modification_type"

  values_color_arrow <- c("CATALYSIS")

  colores_arrow <-c("#848484")

  setEdgeTargetArrowColorMapping(column_line_color_arrow,values_color_arrow,colores_arrow,
                                 mapping.type = "d")

  ######COLORES DE LAS LINEAS
  ###Mal, tiene que hacerse uno por uno como
  ### en la arrow T para inhibition

  setEdgeColorDefault("#848484")
  #column_line_color_edge <- "modification_type"
  #values_line_color_edge <- c("PHYSICAL_STIMULATION","INHIBITION")

  #line_color_edge <- c("#48c4dc","#d80c0c")


  #setEdgeColorMapping(column_line_color_edge,values_line_color_edge,line_color_edge,
  #                    mapping.type = "d")


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

  ###Remover todas las flechas que apunten a los auxiliares

  edges_auxiliares <- edgedata[grepl("Aux*", edgedata$target),"name"]
  setEdgeTargetArrowShapeBypass(edges_auxiliares,"NONE")


##########Obtener los edges adyacentes a catalisis y hacer que tengan forma
  ###open circle
  edgedata <- getTableColumns("edge")

  edges_catalysis <- filter(edgedata, modification_type == "CATALYSIS")

  ##obtener nodos que son source de catalisis
  nodos_catalisis <- edges_catalysis$source


  ###Obtener las edges que tienen como target a los nodos_catalisis

  edges_target_catalisis <- filter(edgedata, target %in% nodos_catalisis)

  ###Obtener solo los auxiliares que no tengan "CATALYSIS"

  auxiliares_edges_target_catalsis <- edges_target_catalisis[!grepl("CATALYSIS", edges_target_catalisis$modification_type),]

  nombres_edges_target_catalisis <- as.character(auxiliares_edges_target_catalsis$name)

  setEdgeTargetArrowShapeBypass(nombres_edges_target_catalisis,"OPEN_CIRCLE")

##############

##########Convertir las lineas de reacciones de transcripción a tipo “marquee dash dot”

  edges_transcription <- edgedata[edgedata$`shared interaction` == "TRANSCRIPTION","name"]

  setEdgeLineStyleBypass(edges_transcription, "MARQUEE_DASH_DOT" )


  #############


  ####Convertir las líneas de modificacion que inciden sobre una reacción de transcripcion
  ####a tipo “Contiguous Arrow”. Solo las que inciden sobre reacciones de transcripcion,
  ####las otras dejarlas como están.

  ######Obtener nodos auxiliares adyacentes a transcripcion
  edges_transcription <- edgedata[edgedata$`shared interaction` == "TRANSCRIPTION",]

  nodos_aux_transcripcion <- edges_transcription[grepl("Aux_", edges_transcription$target),"target"]
  ####Obtener edges que tienen como target a este nodo y que no sean los primeros edges
  edges_transcription_target <- edgedata[edgedata$target %in% nodos_aux_transcripcion,]

  nombres_edges_transcription_original <- as.character(edges_transcription$name)
  edges_transcription_target <- edges_transcription_target[!edges_transcription_target$name %in% nombres_edges_transcription_original,"name" ]

  setEdgeLineStyleBypass(edges_transcription_target, "CONTIGUOUS_ARROW" )

##############

  #Terminar en punta (Target arrow shape) tipo “T”
  #todas las modificaciones de tipo “INHIBITION”.
  #Cuidado de mantener el color de la línea en la punta.



  edges_inhibition <- filter(edgedata, modification_type == "INHIBITION")

  #Verificar si hay inhibiciones

  if (length(rownames(edges_inhibition))>=1){



  nodos_originales_inhibition <- filter(modifications, modification_type =="INHIBITION")
  nodos_originales_inhibition <- unique(as.character(nodos_originales_inhibition$Object))
  #Obtener las reacciones re_Re
  edges_inhibition <- filter(edges_inhibition, reaction_id_direccion == "re_Re")

  ###Obtener los nodos auxuliares de las reacciones re_Re
  nodos_aux_inhibition <- unique(as.character(edges_inhibition$target))

  ###Obtener las edges que tengan como target los nodos aux inhibition y como source
  ###los nodos originales inhibition

  edges_inhibition_1 <- filter(edgedata, target %in% nodos_aux_inhibition &
                                 source %in% nodos_originales_inhibition)


  edges_inhibition_1 <- as.character(edges_inhibition_1$name)
  setEdgeTargetArrowShapeBypass(edges_inhibition_1,"T")
  setEdgeColorBypass(edges_inhibition_1,"#d80c0c")
  setEdgeTargetArrowColorBypass(edges_inhibition_1,"#d80c0c")
  } else{
    print ("No hay Inhibiciones")
  }


  #########Physical Stimulations ####MAL, SE COMENATARA PARA BUSCAR EDGES QUE
  ########TAMBIEN SEAN INHIBIDOS Y COLOREARLOS DE AMARILLO

  edges_stimulation <- filter(edgedata, modification_type == "PHYSICAL_STIMULATION")

  #Verificar si hay inhibiciones

  if (length(rownames(edges_stimulation))>=1){

    nodos_originales_stimulation <- filter(modifications, modification_type =="PHYSICAL_STIMULATION")
    nodos_originales_stimulation <- unique(as.character(nodos_originales_stimulation$Object))
    #Obtener las reacciones re_Re
    edges_stimulation <- filter(edges_stimulation, reaction_id_direccion == "re_Re")

    ###Obtener los nodos auxuliares de las reacciones re_Re
    nodos_aux_stimulation <- unique(as.character(edges_stimulation$target))

    ###Obtener las edges que tengan como target los nodos aux stimulation y como source
    ###los nodos originales stimulation

    edges_stimulation_1 <- filter(edgedata, target %in% nodos_aux_stimulation &
                                   source %in% nodos_originales_stimulation)


    edges_stimulation_1 <- as.character(edges_stimulation_1$name)
    setEdgeTargetArrowShapeBypass(edges_stimulation_1,"DELTA")
    setEdgeColorBypass(edges_stimulation_1,"#0099cc")
    setEdgeTargetArrowColorBypass(edges_stimulation_1,"#0099cc")
  } else{
    print ("No hay modificaciones Physical Stimulation")
  }

  #####################

  ######Terminar en punta (Target arrow shape) tipo “Cross Open Delta”
  #####todas las líneas que apunten a productos de reacciones de tipo “TRANSPORT”.

  edges_transport <- filter(edgedata, `shared interaction` == "TRANSPORT")
  if (length(rownames(edges_transport))>=1){
    ###obtener las que van en direccion aux_reactivo
    edges_transport <- filter(edges_transport, direccion == "aux_reactivo")
    edges_transport <- as.character(edges_transport$name)
    setEdgeTargetArrowShapeBypass(edges_transport,"CROSS_OPEN_DELTA")



  } else{

  }



  ####PHYSYCAL STIMULATION BIEN( colorea las reacciones)
  #edges_stimulation <- filter(edgedata, modification_type == "PHYSICAL_STIMULATION")
  #Verificar si hay STIMULATIONS


  #if (length(rownames(edges_stimulation))>=1){

    #nodos_originales_stimulation <- edges_stimulation

    #nodos_originales_stimulation <- unique(as.character(nodos_originales_stimulation$Object))
    #Obtener las reacciones re_Re
    #edges_stimulation <- filter(edges_stimulation, reaction_id_direccion == "re_Re")

    ###Obtener los nodos auxuliares de las reacciones re_Re
    #nodos_aux_stimulation <- unique(as.character(edges_stimulation$target))

    ###Obtener las edges que tengan como target los nodos aux stimulation y como source
    ###los nodos originales stimulation

    #edges_stimulation_1 <- filter(edgedata, target %in% nodos_aux_stimulation &
                                    #source %in% nodos_originales_stimulation)


    #edges_stimulation_1 <- as.character(nodos_originales_stimulation$name)
    #setEdgeTargetArrowShapeBypass(edges_stimulation_1,"DELTA")
    #setEdgeColorBypass(edges_stimulation_1,"#0099cc")
    #setEdgeTargetArrowColorBypass(edges_stimulation_1,"#0099cc")

    ###Remover todas las flechas que apunten a los auxiliares

    #edges_auxiliares_stimulation <- edges_stimulation[grepl("Aux*", edges_stimulation$target),"name"]
    #setEdgeTargetArrowShapeBypass(edges_auxiliares_stimulation,"NONE")

  #} else{
  #  print ("No hay modificaciones Physical Stimulation")
  #}


  ##########VERIFICAR AQUELLAS MODIFICACIONES QUE SEAN
  ##########INHIBITION Y PHYSYCAL STIMULATION EN LA MISMA REACCION

  modifications_mod <- read.delim("modification.txt", header = FALSE,
                              col.names = c("reaction","modification_type"
                                            ,"Object","Reference"), stringsAsFactors=FALSE)

  modificaciones_inhibition <- filter(modifications_mod, modification_type =="INHIBITION")
  modificaciones_stimulation <- filter(modifications_mod, modification_type =="PHYSICAL_STIMULATION")

  modificaciones_dobles <- inner_join(modificaciones_stimulation, modificaciones_inhibition, by="reaction" )

  modificaciones_dobles_nombres <- as.character(unique(modificaciones_dobles$reaction))


  if ((length(modificaciones_dobles_nombres))>=1){

    modificaciones_dobles_RE <- paste(modificaciones_dobles_nombres, "Re", sep="_")
    modificaciones_dobles_PR <- paste(modificaciones_dobles_nombres, "Pr", sep="_")

    modificaciones_dobles_reaction <- c(modificaciones_dobles_PR, modificaciones_dobles_RE)

    modificaciones_dobles <- filter(edgedata, reaction_id_ %in% modificaciones_dobles_reaction)
    modificaciones_dobles_object <- as.character(modificaciones_dobles$reaction_id_)


    #####Obtener los nodos originales
    nodos_originales_mod_dobles <- filter(modifications, reaction %in% modificaciones_dobles_object)
    nodos_originales_mod_dobles<- unique(as.character(nodos_originales_mod_dobles$Object))

    ####Obtener los nodos auxiliares que sean target de estas reacciones Re

    edges_aux_mod_dobles <- filter(modificaciones_dobles, reaction_id_direccion == "re_Re")
    ###Obtener los nodos auxuliares de las reacciones re_Re
    nodos_aux_mod_dobles <- unique(as.character(edges_aux_mod_dobles$target))

    ###Obtener las edges que tengan como target los nodos_aux_mod_dobles y como source
    ###los nodos_originales_mod_dobles

    edges_mod_dobles_1 <- filter(edgedata, target %in% nodos_aux_mod_dobles &
                                   source %in% nodos_originales_mod_dobles)

    edges_mod_dobles_1 <- as.character(edges_mod_dobles_1$name)





    setEdgeColorBypass(edges_mod_dobles_1,"#ffcc33")
    setEdgeTargetArrowColorBypass(edges_mod_dobles_1,"#ffcc33")


  }


  ###Remover todas las flechas salgan de los auxiliares con _Pr





  ################

  # if(!(missing(html))){
  #  if(html){
  #suppressWarnings(RegulonDBCytoscape::exportToHTML())
  # }
  #}


}
