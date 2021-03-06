% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.assay.file.R
\name{read.assay.file}
\alias{read.assay.file}
\title{read the assay file under the dir}
\usage{
read.assay.file(
  dir = "./",
  pattern,
  plate.num = 4,
  skip = 11,
  plate.pattern = "plate",
  out.plate = F,
  snp.name.ix = 3,
  plate.ix = 4,
  gene.name.ix = 2
)
}
\arguments{
\item{dir}{all raw assay file are under this dir.}

\item{pattern}{the prefix or extension name of the all raw files.}

\item{plate.num}{the number of plates used in the assay.}

\item{skip}{the skip line in the head of each assay file.}

\item{plate.pattern}{the prefix name of the plate (default:"plate").}

\item{out.plate}{logical. if \code{TRUE}, the output will add a col contain the plate.num of each sample.}

\item{snp.name.ix}{the index of snp name after this function split the file name (default:3).}

\item{plate.ix}{the index of plate number after this function split the file name (default:4).}

\item{gene.name.ix}{the index of gene name after this function split the file name (default:2).}
}
\value{
a data frame with the col for sample and the row for snp.loci:
\item{output}{function output}
}
\description{
\code{read.assay.file} read the assay files under the dir, and the output can be the input of function \code{from.assay.to.ped}.
}
\details{
This function read the raw assay files and merge them. The output of this function is 
a data.frame. The each row of it is a sample. And each col is a snp loci.
}
\examples{
data.total.1536 <- read.assay.file(dir = "./", pattern = "AZ")



}
\author{
Liangjie Liu <liuliangjie@sjtu.edu.cn>
}
