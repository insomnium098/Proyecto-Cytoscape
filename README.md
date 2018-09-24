# CytoscapeRegulon

This R package creates a Cytoscape network from GENSOR unit files from Regulon DB. 

Each transcription factor folder needs to have 5 tab delimited files:
complexes.txt, modification.txt, objects.txt, reactants_products.txt and reactions.txt. 

Also, Cytoscape needs to be open at the same time as the package is running.


-----INSTALLATION

A) Directly from github

library(devtools)


install_github("insomnium098/Proyecto-Cytoscape/RegulonDBCytoscape")

OR


B) Using binary package (only for mac OS)

Download the binary file from https://github.com/insomnium098/Proyecto-Cytoscape/releases

In RStudio select Tools > Install packages > And in package archive select the downloaded package



----- EXAMPLE USAGE:

x <- c("/Users/daniel/Proyecto-Cytoscape/SRsGUs/LacI")

cytoscapeRegulon(x)
