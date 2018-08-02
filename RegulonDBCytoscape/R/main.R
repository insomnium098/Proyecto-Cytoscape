library(dplyr)

#x es la ruta de la carpeta con los archivos
x <- c("AdiY")

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


  ####



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


  #Llamar RCy3
  createNetworkFromDataFrames(df_nodes,df_edges, title=x,
                              collection=x)


}
