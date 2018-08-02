library(dplyr)
library(RCy3)
library(plyr)

#x es la ruta de la carpeta con los archivos
x <- c("LacI")

cytoscapeRegulon <- function(x){
  #x es el folder con los archivos
  setwd(x)
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





  ####Complejos

  #Paso1: obtener los complejos

  complejos <- unique(complexes$id)

  #Paso2: extraer cada id de complejo

  for (i in complejos){
    complejo_df <- complexes[complexes$id==i,]
    nombre_complejo <- i
    #crear data frame para cada integrante del complejo
    for (j in seq_along(1:length(rownames(complejo_df)))){
      integrante_complejo <- complejo_df[j,]
      nombre_integrante_complejo <- as.character(integrante_complejo$type)
      df_integrante_complejo <- data.frame("id" = paste0(
        nombre_integrante_complejo," "),
        "reactant_type" = "COMPLEX_NODE",
        "name" = i)
      #HACER UN DF_NODES CON COMPLEJOS
      if (!exists("df_nodes_complex")){
        df_nodes_complex <- df_integrante_complejo
      }


      #Si el dataset existe, se une
      if (exists("df_nodes_complex")){
        temp_dataset <- df_integrante_complejo
        #dataset<-rbind(dataset, temp_dataset)
        df_nodes_complex <- rbind(df_nodes_complex,temp_dataset)
        rm(temp_dataset)
      }

    }

  }

  #El primer row se repite, se eliminara
  df_nodes_complex <- df_nodes_complex[2:length(rownames(df_nodes_complex)),]
  #Unir el df_nodes_complex con el df_nodes

  df_nodes <- rbind(df_nodes, df_nodes_complex)


  ######Modificar el df de los complejos con los nuevos nombres con espacio

  nombres_complejos <- as.character(complexes$type)
  nombres_complejos <- paste0(nombres_complejos, ' ')

  complexes$type <- nombres_complejos

  ###Cambiar los colnames de complejes para posteriormente hacerlos edges
  colnames(complexes) <- c("reaction_id_","target","source")

  complexes$interaction <- rep("COMPLEX_INTERACTION",
                               length(rownames(complexes)))

  ##unir complexes con df_edges

  df_edges <- rbind.fill(df_edges,complexes)




  #Llamar RCy3
  createNetworkFromDataFrames(df_nodes,df_edges, title=x,
                              collection=x)



  ####Agrupar moleculas de los complejos

  nodedata <- getTableColumns("node")
  nodes_complexes <- nodedata[grep("COMPLEX_NODE", nodedata$reactant_type), ]


  #####

  # Seleccionar los complejos, uno por uno
  for (g in unique(as.character(complexes$reaction_id_))){



    deltacatnodes <- df_nodes[grep(g, df_nodes$name), "id"]
    ##OBTENER EL SUID
    SUIDS <- filter(nodedata,`shared name` %in% deltacatnodes)
    SUIDS <- as.character(SUIDS$SUID)
    ####
    nodos_suids <- selectNodes(SUIDS, preserve=FALSE)
    ###Seleccionar el primer nodo vecino
    primer_vecino <- selectFirstNeighbors()
    primer_vecino <- setdiff(primer_vecino,nodos_suids)
    primer_vecino <- filter(nodedata, SUID == primer_vecino)
    primer_vecino <- as.character(primer_vecino$`shared name`)

    ##Crear grupo del complejo
    createGroup(primer_vecino)
    collapseGroup(primer_vecino)
  }
  ########

  #######Forma de los nodos
  #######
  column <- 'reactant_type'
  values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
               'RNA',"GENE","COMPLEX","PROTEIN_TF","COMPLEX_NODE")
  shapes <- c ('ELLIPSE', 'ROUND_RECTANGLE',
               'PARALLELOGRAM',"RECTANGLE","OCTAGON",
               "ROUND_RECTANGLE","ELLIPSE")
  setNodeShapeMapping (column, values, shapes)

  #####Color de los nodos
  column <- 'reactant_type'
  values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
               'RNA',"GENE","COMPLEX","PROTEIN_TF","COMPLEX_NODE")
  colors <- c("#b6bd7b","#b6bd7b","#ffbc00","#ffbc00","#b6bd7b","#4881a6",
              "#b6bd7b")
  setNodeColorMapping (column, values, colors,
                       mapping.type = "d")
  ###Arrow shapes

  #getArrowShapes()
  #?setEdgeLineStyleMapping

  column_arrow <- "interaction"
  values_arrow <- c("STATE_TRANSITION","TRANSCRIPTION","TRANSLATION",
                    "TRANSPORT","COMPLEX_INTERACTION")

}
