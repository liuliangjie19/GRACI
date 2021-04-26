#' @title the plot function for the Object of class GGC.data
#'
#' @description
#' \code{plot} For the object of class GGC.data, plot the object as a network.
#'
#' @details
#' This function use the output of the function \code{gene.gene.causal}. 
#' The output image of this function is a network. Every node in the network
#' is a gene, and every edge in the network is a causal relationship between 
#' two genes. And the node can be plot in different color or shape.
#' 
#' @param ggc output of the function \code{gene.gene.causal}.
#' @param color.method (optional) a character string matching one of the following items: snp and out.
#' @param color.threshold the threshold to define the color of each node in the network.
#' @param shape.method (optional) a character string matching one of the following items: snp and out.
#' @param shape.threshold the threshold to define the shape of each node in the network.
#' @method plot GGC.data
#' @export
#' @return a network plot contain all the gene and the relationship between them.
#' @author Liangjie Liu <liuliangjie@@sjtu.edu.cn>
#' @examples
#' plot(temp.ggc, color.method = "snp", shape.method = "out")
#' 
#' 

plot.GGC.data <- function(ggc, color.method = c("snp", "out"), color.threshold = 5, shape.method = c("snp", "out"), shape.threshold = 5){
  #require(igraph)
  color.method <- match.arg(color.method)
  shape.method <- match.arg(shape.method)

  graph.data <- ggc$gg.causal
  graph.vertices <- ggc$gene.table
  ngenes <- dim(graph.vertices)[1]

  graph.vertices$out.edge <- rep(0, ngenes)
  for (gene in graph.vertices[,1]){
    graph.vertices[which(graph.vertices[,1]==gene),]$out.edge <- sum(graph.data[,1]==gene)
  }

  graph.vertices$color <- rep("", ngenes)
  if (color.method=="snp"){
    #graph.vertices[which(graph.vertices$Snps.Freq>color.threshold),]$color <- "red"
    #graph.vertices[which(graph.vertices$Snps.Freq<=color.threshold),]$color <- "blue"
    graph.vertices$color <- ifelse(graph.vertices$Snps.Freq>color.threshold, "red", "yellow")
  }
  else if (color.method=="out"){
    #graph.vertices[which(graph.vertices$out.edge>color.threshold),]$color <- "red"
    #graph.vertices[which(graph.vertices$out.edge<=color.threshold),]$color <- "blue"
    graph.vertices$color <- ifelse(graph.vertices$out.edge>color.threshold, "red", "yellow")
  }

  graph.vertices$shape <- rep("", ngenes)
  if (shape.method=="snp"){
    #graph.vertices[which(graph.vertices$Snps.Freq>shape.threshold), ]$shape <- "rectangle"
    #graph.vertices[which(graph.vertices$Snps.Freq<=shape.threshold),]$shape <- "pie"
    graph.vertices$shape <- ifelse(graph.vertices$Snps.Freq>shape.threshold, "rectangle", "circle")
  }
  else if (shape.method=="out"){
    #graph.vertices[which(graph.vertices$out.edge>shape.threshold), ]$shape <- "rectangle"
    #graph.vertices[which(graph.vertices$out.edge<=shape.threshold),]$shape <- "pie"
    graph.vertices$shape <- ifelse(graph.vertices$out.edge>shape.threshold, "rectangle", "circle")
  }

  print(graph.vertices)
  print(graph.data)

  network.ggc <- igraph::graph_from_data_frame(graph.data, directed = T, vertices = graph.vertices)
  plot(network.ggc, vertex.shapes = V(network.ggc)$shape, vertex.color = V(network.ggc)$color)
}
