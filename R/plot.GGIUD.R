#' @title the plot function for the Object of class GGIUD
#'
#' @description
#' \code{plot} For the object of class GGIUD, plot the object as a network.
#'
#' @details
#' This function use the output of the function \code{gene.gene.interaction}. 
#' The output image of this function is a network. Every node in the network
#' is a gene, and every edge in the network is a causal relationship between 
#' two genes. And the node can be plot in different color or shape.
#' 
#' @param ggiud output of the function \code{gene.gene.interaction}.
#' @param p.threshold the threshold for p.value, interactions with p.value less than p.threshold are vaild.
#' @param color.method (optional) a character string matching one of the following items: snp and out degree.
#' @param color.threshold the threshold to define the color of each node in the network.
#' @param shape.method (optional) a character string matching one of the following items: snp and out degree.
#' @param shape.threshold the threshold to define the shape of each node in the network.
#' @method plot GGIUD
#' @export
#' @return a network plot contain all the gene and the relationship between them.
#' @author Liangjie Liu <liuliangjie@@sjtu.edu.cn>
#' @examples
#' plot(temp.ggiud, color.method = "snp", shape.method = "out")
#' 

plot.GGIUD <- function(ggiud, p.threshold = 0.05, color.method = c("snp", "degree"), color.threshold = 5, shape.method = c("snp", "degree"), shape.threshold = 5) {
  require(igraph)

  color.method <- match.arg(color.method)
  shape.method <- match.arg(shape.method)

  graph.data <- ggiud$GGI.ud[ggiud$GGI.ud$P.val<p.threshold, ]
  graph.vertices <- ggiud$GENE.table
  ngenes <- dim(graph.vertices)[1]

  graph.vertices$degree <- rep(0, ngenes)
  for (gene in graph.vertices[,1]){
    graph.vertices[which(graph.vertices[,1]==gene),]$degree <- sum(graph.data[,1]==gene) + sum(graph.data[,2]==gene)
  }
  graph.vertices$color <- rep("", ngenes)
  if (color.method=="snp"){
    graph.vertices$color <- ifelse(graph.vertices$Snps.Freq > color.threshold, "red", "yellow")
  }
  else if (color.method=="degree"){
    graph.vertices$color <- ifelse(graph.vertices$degree > color.threshold, "red", "yellow")
  }

  graph.vertices$shape <- rep("", ngenes)
  if (shape.method=="snp"){
    graph.vertices$shape <- ifelse(graph.vertices$Snps.Freq>shape.threshold, "rectangle", "circle")
  }
  else if (shape.method=="degree"){
    graph.vertices$shape <- ifelse(graph.vertices$degree > shape.threshold, "rectangle", "circle")
  }

  print(graph.vertices)
  print(graph.data)

  network.ggiud <- igraph::graph_from_data_frame(graph.data, directed = F, vertices = graph.vertices)
  plot(network.ggiud, vertex.shapes = V(network.ggiud)$shape, vertex.color = V(network.ggiud)$color)
}
