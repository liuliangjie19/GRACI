% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.GGIUD.R
\name{plot.GGIUD}
\alias{plot.GGIUD}
\title{the plot function for the Object of class GGIUD}
\usage{
\method{plot}{GGIUD}(
  ggiud,
  p.threshold = 0.05,
  color.method = c("snp", "degree"),
  color.threshold = 5,
  shape.method = c("snp", "degree"),
  shape.threshold = 5
)
}
\arguments{
\item{ggiud}{output of the function \code{gene.gene.interaction}.}

\item{p.threshold}{the threshold for p.value, interactions with p.value less than p.threshold are vaild.}

\item{color.method}{(optional) a character string matching one of the following items: snp and out degree.}

\item{color.threshold}{the threshold to define the color of each node in the network.}

\item{shape.method}{(optional) a character string matching one of the following items: snp and out degree.}

\item{shape.threshold}{the threshold to define the shape of each node in the network.}
}
\value{
a network plot contain all the gene and the relationship between them.
}
\description{
\code{plot} For the object of class GGIUD, plot the object as a network.
}
\details{
This function use the output of the function \code{gene.gene.interaction}. 
The output image of this function is a network. Every node in the network
is a gene, and every edge in the network is a causal relationship between 
two genes. And the node can be plot in different color or shape.
}
\examples{
plot(temp.ggiud, color.method = "snp", shape.method = "out")

}
\author{
Liangjie Liu <liuliangjie@sjtu.edu.cn>
}
