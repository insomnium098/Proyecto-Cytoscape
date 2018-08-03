% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{cytoscapeRegulon}
\alias{cytoscapeRegulon}
\title{Plot MDS for a DESeq2 object}
\usage{
cytoscapeRegulon(x)
}
\arguments{
\item{x}{A directory with complexes.txt, modification.txt, objects.txt, reactants_products.txt, reactions.txt}

\item{html}{Boolean if True it generates an HTML file.}
}
\description{
This function takes all the output files from "blabla" and conectes them with help of Cytoscape
}
\examples{
##No RUN##
setwd("myDirectoryofTFs/")
x<-LacI
cytoscapeRegulon(x, html=TRUE)

}
\author{
Antonio Martinez, Luis Pedro Iniguez, Servando Ramirez, Tobias Portillo
}