library(dplyr)


########## Procesamiento de reactants_products.txt y reactions.txt
##########  El producto es df_edges

reactants_products <- read.delim("reactants_products.txt", header = FALSE,
                                 col.names = c("id","name","type") )


reactions <- read.delim("reactions.txt", header = FALSE,
                        col.names = c("id","reaction_type"
                                      ,"direction",
                                      "ecocyc_name_id"))


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

#Eliminar dataframes que ya no se usan

rm(df_product, df_reactant, df_reactant_product, reactants_products, reactions)

#####
#####
