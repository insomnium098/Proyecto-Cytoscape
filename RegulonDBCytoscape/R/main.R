library(dplyr)
library(RCy3)

#x es la ruta de la carpeta con los archivos
x <- c("AtoC")

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



  #####

  ###Separar matriz reactant_product por producto y reaccion para
  ## hacer una tercer data.frame formateada para RCY3

  df_reactant <- subset(reactants_products, type == "reactant")
  df_product <- subset(reactants_products, type == "product")
  df_reactant_product <- merge(df_reactant, df_product, by="id")
  df_reactant_product <- df_reactant_product[,c(1,2,4)]
  colnames(df_reactant_product) <- c("id","reactant_name",
                                     "product_name")

  ###Añadir datos de reactions en dataframe con edges final
  df_edges <- left_join(df_reactant_product, reactions, by = "id")

  colnames(df_edges) <- c("reaction_id","source", "target","interaction",
                          "direction","ecocyc_name")
  ##
  df_edges$reaction_id_ <- df_edges$reaction_id
  df_edges <- df_edges[,2:length(colnames(df_edges))]
  ##
  #Eliminar dataframes que ya no se usan

  rm(df_product, df_reactant, df_reactant_product, reactants_products, reactions)


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
        #Se añadira un espacio para que el nodo nuevo no tenga interaccion
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



  #Llamar RCy3
  createNetworkFromDataFrames(df_nodes,df_edges, title=x,
                              collection=x)

  #######Forma de los nodos
  #######
  column <- 'reactant_type'
  values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
               'RNA',"GENE","COMPLEX","PROTEIN_TF")
  shapes <- c ('ELLIPSE', 'ROUND_RECTANGLE',
               'PARALLELOGRAM',"RECTANGLE","OCTAGON",
               "ROUND_RECTANGLE")
  setNodeShapeMapping (column, values, shapes)

  #####Color de los nodos
  column <- 'reactant_type'
  values <- c ('SIMPLE_MOLECULE',  'PROTEIN',
               'RNA',"GENE","COMPLEX","PROTEIN_TF")
  colors <- c("#b6bd7b","#b6bd7b","#ffbc00","#ffbc00","#b6bd7b","#4881a6")
  setNodeColorMapping (column, values, colors,
                       mapping.type = "d")









}
