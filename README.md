# CytoscapeRegulon

This R package creates a Cytoscape network from GENSOR unit files from Regulon DB. 

Each transcription factor folder needs to have 5 tab delimited files:
complexes.txt, modification.txt, objects.txt, reactants_products.txt and reactions.txt. 

Also, Cytoscape needs to be open at the same time as the package is running.

Example usage:

x <- c("/Users/daniel/Proyecto-Cytoscape/SRsGUs/LacI")

cytoscapeRegulon(x)
