#' Plot MDS for a DESeq2 object
#'
#' @description This function takes all the output files from "blabla" and conectes them with help of Cytoscape
#' @author Antonio Martinez, Luis Pedro Iniguez, Servando Ramirez, Tobias Portillo
#' @param x A directory with complexes.txt, modification.txt, objects.txt, reactants_products.txt, reactions.txt
#' @param html Boolean if True it generates an HTML file.
#' @export
#' @import RCy3


generateImage <- function(){
  full.path=paste(getwd(),'graph_image',sep='/')
  exportImage(filename=full.path, type = 'JPG')
}

copytToWd <- function(){
  if (dir.exists(paste(getwd(),'/websiteRegulon'))){
    TRUE
  }else{
    file.copy(paste(find.package("RegulonDBCytoscape"),'websiteRegulon',sep='/'),paste(getwd()), recursive=TRUE)
    TRUE
  }
}

createStyle <- function(){
  RCy3::exportVisualStyles(paste(getwd(),'style',sep='/'),type="JSON")
  firstLine="var styles ="
  fConn = file(paste(getwd(),'style.json',sep='/'), 'r+')
  lines <- readLines(fConn)
  text <- c(firstLine,lines)
  write(text,file=paste(getwd(),'websiteRegulon','data','styles.js',sep='/'),append=FALSE)
}

createJavascript <- function(){
  RCy3::exportNetwork(paste(getwd(),'data',sep='/'),type="CYJS")
  line="var networks = {'lacI_testnetwrk.txt':"
  line2="}"
  fConn = file(paste(getwd(),'data.cyjs',sep='/'), 'r+')
  Lines <- readLines(fConn)
  text <- c(line,Lines,line2)
  write(text,file=paste(getwd(),'websiteRegulon','data','networks.js',sep='/'),append=FALSE)
}

exportToHTML <- function(){
  if (copytToWd()){
    createJavascript()
    createStyle()
    generateImage()
  }
}

