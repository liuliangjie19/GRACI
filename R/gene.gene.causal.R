#' @title calculate the causal relationship between the interacting gene pairs
#'
#' @description
#' \code{gene.gene.causal} calculate the causal relationship between the interacting gene pairs, although you can calculate the causal relationship between any gene pairs, that doesn't make senses.
#'
#' @details
#' After the function \code{gene.gene.interaction} calculates the 
#' interaction between genes, the causal relationship between gene 
#' pairs with obvious interaction is inferred. When f<0, gene 1 leads 
#' to gene 2, and f>0 is the opposite.
#' 
#' @param data.gguid the ouput of the function \code{gene.gene.interaction}, a list of the class GGIUD.
#' @param n.times the times of the permutation.
#' @param sample.size the sample.size of each time of permutation
#' @param ref (optional) a character string matching one of the following items: Gaussian and uniform. the data distribution of your data in function \code{igci}
#' @param est (optional) a character string matching one of the following items: entropy and intergral. the estimator used in function \code{igci}
#' @param interaction.p the p.value threshold of the interacting relationship, default: 0.05.
#' @param fex logical. if \code{TRUE}, using the function \code{fex} to reduce the dimensions of SnpMatrix.
#' @export
#' @return a GGC.data list contain 6 item:
#' \item{SnpMatrix}{function output}
#' \item{gene.info}{gene info}
#' \item{sample.info}{sample info}
#' \item{gg.causal}{gene gene causal table}
#' \item{gene.table}{gene table}
#' \item{gg.f}{gene f.value table}
#' @author Liangjie Liu <liuliangjie@@sjtu.edu.cn>
#' @examples
#' temp.ggc <- gene.gene.causal(temp.ggiud, n.times = 1000, sample.size = 700, fex = TRUE)
#' 

gene.gene.causal <- function(data.ggiud, n.times = 1000, sample.size = 100,
                             ref = c("Gaussian", "uniform"), est = c("entropy", "intergral"),
                             interaction.p = 0.05, fex = T, pos.use = T) {
  #if (!inherits(data.ggigci, "GGIGCI.data")){
  #  stop("data must be one of the class 'GGIGCI.data'(the output of the function impute.res).")
  #}
  if (!inherits(data.ggiud, "GGIUD")){
    stop("data must be one of the class 'GGIUD'(the output of the function gene.gene.interaction).")
  }

  ref <- match.arg(ref)
  est <- match.arg(est)
  #print(paste0(ref, "xiix"))
  #print(est)
  gene.table <- data.ggiud$GENE.table
  snp.info <- data.ggiud$SnpMatrix
  gene.info <- data.ggiud$gene.info
  sample.info <- data.ggiud$sample.info
  phen <- as(sample.info[,6], "numeric")
  selected.gene.pairs <- data.ggiud$GGI.ud[which(data.ggiud$GGI.ud$P.val<=interaction.p), ]
  gene.f.table <- data.frame(gene1 = selected.gene.pairs[, 1], gene2 = selected.gene.pairs[, 2],
                             f1 = rep(NA, length(selected.gene.pairs[,1])), f2 = rep(NA, length(selected.gene.pairs[, 2])),
                             P.value = rep(0, length(selected.gene.pairs[, 1])), causal.direction = rep(0, length(selected.gene.pairs[, 1])))
  gene.causal.table <- data.frame(from = selected.gene.pairs[, 1], to = selected.gene.pairs[, 2],
                                  F.value = rep(0, length(selected.gene.pairs[, 1])))
  n.pairs <- dim(selected.gene.pairs)[1]
  for(i in seq_len(n.pairs)){
    #print(selected.gene.pairs[i,])
    gene1 <- selected.gene.pairs[i, 1]
    gene2 <- selected.gene.pairs[i, 2]
    snp.gene1 <- gene.info[which(gene.info[,2]==gene1), 3]
    snp.gene2 <- gene.info[which(gene.info[,2]==gene2), 3]
    snp.info.gene1 <- as(snp.info[, colnames(snp.info)%in%snp.gene1], "numeric")
    snp.info.gene2 <- as(snp.info[, colnames(snp.info)%in%snp.gene2], "numeric")

    if (fex){
      #x.gene1 <- fex(cbind(phen, snp.info.gene1))
      #x.gene2 <- fex(cbind(phen, snp.info.gene2))
      if(pos.use){
        x.gene1 <- fex(snp.info.gene1, gene.info[which(gene.info[,2]==gene1), ])
        x.gene2 <- fex(snp.info.gene2, gene.info[which(gene.info[,2]==gene2), ])
      }
      else {
        x.gene1 <- fex(snp.info.gene1)
        x.gene2 <- fex(snp.info.gene2)
      }

      F.matrix <- diag(0, nrow = dim(x.gene1)[2], ncol = dim(x.gene2)[2])
      rownames(F.matrix) <- colnames(x.gene1)
      colnames(F.matrix) <- colnames(x.gene2)
      for(n in seq_len(dim(F.matrix)[1])) {
        for(m in seq_len(dim(F.matrix)[2])) {
          pairs.data <- data.frame(gene1 = x.gene1[, n], gene2 = x.gene2[, m])
          #fex.temp.result <- IGCIfor2value(pairs.data, times = n.times, sample.size = sample.size, refMeasure = ref, estimator = est)
          F.matrix[n, m] <- igci(pairs.data[,1], pairs.data[, 2], refMeasure = ref, estimator = est)
        }
      }
      F1 <- sum(F.matrix<0)
      F2 <- sum(F.matrix>0)
      gene.f.table[i, 3:5] <- c(F1, F2, NA)
      if (F1>F2){
        gene.f.table[i, 5] <- F2/(F1+F2)
        gene.f.table[i, 6] <- "gene1->gene2"
      }
      else if (F1<F2) {
        gene.f.table[i, 5] <- F1/(F1+F2)
        gene.f.table[i, 6] <- "gene2->gene1"
      }
    }
    else {
      x.gene1 <- apply(cbind(phen, 2*base3to10(snp.info.gene1)), 1, sum)
      x.gene2 <- apply(cbind(phen, 2*base3to10(snp.info.gene2)), 1, sum)
      pairs.data <- data.frame(gene1 = x.gene1, gene2 = x.gene2)
      temp.result <- IGCIfor2value(pairs.data, times = n.times, sample.size = sample.size, refMeasure = ref, estimator = est)

      gene.f.table[i, 3:6] <- c(temp.result$f.mean[1], temp.result$f.mean[2], temp.result$p.value, temp.result$causal.direction)
    }
    if (gene.f.table[i, 6]=="gene1->gene2"){
      gene.causal.table[i, ] <- c(gene.f.table[i, 1], gene.f.table[i, 2], gene.f.table[i, 3])
    }
    else if (gene.f.table[i, 6]=="gene2->gene1") {
      gene.causal.table[i, ] <- c(gene.f.table[i, 2], gene.f.table[i, 1], gene.f.table[i, 4])
    }
  }
  out <- list(SnpMatrix = snp.info, gene.info = gene.info, sample.info = sample.info,
              gg.causal = gene.causal.table, gene.table = gene.table, gg.f = gene.f.table)
  if (fex){
    out$method <- "Fexpand"
  }
  else {
    out$method <- "base3to10"
  }

  class(out) <- c("GGC.data", "list")
  return(out)
}

#fex <- function(gene, pos){
#
#  return(apply(gene, 1, sum))
#}
fex <- function(gene, pos=NULL) {
  if (is.null(pos)) {
    pos.v <- (0:(dim(gene)[2]-1))/(dim(gene)[2]-1)
    names(pos.v) <- colnames(gene)
    gene.name <- "gene"
  }
  else {
    pos.v <- as(pos$Position, "numeric")
    names(pos.v) <- pos$SNPnames
    gene.name <- pos$Genenames[1]
  }
  nsnps <- length(pos.v)
  if (nsnps>1){
    idx <- order(pos.v)
    gene <- gene[, idx]
    pos.v <- pos.v[idx]
    pos.v <- (pos.v - pos.v[1])/(pos.v[nsnps]- pos.v[1])
    #print(pos.v)

    eigenval <- prcomp(gene)$sd^2
    eigen.sum <- sum(eigenval)
    tmp <- 0
    n.of.basis <- 0

    for (i in 1:length(eigenval)){
      tmp <- eigenval[i]+tmp
      n.of.basis <- i
      if (tmp >= 0.8*eigen.sum){
        break
      }
    }
    n.of.basis <- floor(n.of.basis/2)*2+1
    #print(n.of.basis)
    frange <- c(pos.v[1], pos.v[length(pos.v)])
    fbasis <- fda::create.fourier.basis(frange, nbasis = n.of.basis)
    phi <- fda::eval.basis(pos.v, fbasis)
    res <- t(MASS::ginv(t(phi)%*%phi)%*%t(phi)%*%t(gene))
  }
  else {
    res <- gene
  }
  col.num <- dim(res)[2]
  #gene.name <- pos$Genenames[1]
  #print(gene.name)
  colnames(res) <- paste(gene.name, seq_len(col.num), sep = ".")

  return(res)
}

base3to10 <- function(gene){
  n.SNP<-ncol(gene)
  long<-nrow(gene)
  tab<-matrix(NA,ncol=n.SNP,nrow=long)
  for (i in seq_len(n.SNP)){
    tab[,i]<-3^(i-1)*gene[,i]
  }
  b10<-apply(tab,1,sum)
  return(b10)
}

# R version of Information Geometry Causal Inference
#                _                           
# platform       x86_64-apple-darwin17.0     
# arch           x86_64                      
# os             darwin17.0                  
# system         x86_64, darwin17.0          
# status                                     
# major          4                           
# minor          0.0                         
# year           2020                        
# month          04                          
# day            24                          
# svn rev        78286                       
# language       R                           
# version.string R version 4.0.0 (2020-04-24)
# nickname       Arbor Day 

# P. Daniusis, D. Janzing, J. Mooij, J. Zscheischler, B. Steudel,
# K. Zhang, B. Scholkopf:  Inferring deterministic causal relations.
# Proceedings of the 26th Annual Conference on Uncertainty in Artificial 
# Intelligence (UAI-2010).  
# http://event.cwi.nl/uai2010/papers/UAI2010_0121.pdf
#######################################################################################
# Inputs:                                                                             #
#   x                       L x 1 observations of x                                   #
#   y                       L x 1 observations of y                                   #
# refMeasure      reference measure to use:                                           #
#   1: uniform                                                                        #
#   2: Gaussian                                                                       #
# estimator       estimator to use:                                                   #
#   1: entropy (eq. (12) in [1]),                                                     #
#   2: integral approximation (eq. (13) in [1]).                                      #
#                                                                                     #
# Outputs:                                                                            #
#   f < 0:          the method prefers the causal direction x -> y                    #
#   f > 0:          the method prefers the causal direction y -> x                    #
#######################################################################################


psi <- function(x) {
  dx <- D(expression(log(gamma(x))), "x")
  return(eval(dx))
}

igci <- function(x, y, refMeasure = "Gaussian", estimator = "entropy") {
  if (refMeasure=="Gaussian" || refMeasure=="g") {
    refMeasure <- 2
  } else if (refMeasure == "uniform" || refMeasure=="u") {
    refMeasure <- 1
  } else {
    stop("refMeasure must be Gaussian or uniform.\n")
  }
  if (estimator=="entropy" || estimator=="e") {
    estimator <- 1
  } else if (estimator == "intergral" || estimator=="i") {
    estimator <- 2
  } else {
    stop("estimator must be extropy or intergral.\n")
  }
  if (!(is.vector(x) & is.vector(y))) {
    stop("x and y must be vector!\n")
  }
  nx <- length(x)
  ny <- length(y)
  if (nx != ny) {
    stop("the length of x and y must be equal!\n")
  } 
  if (nx < 20) {
    stop("Not enough observations in this case!\n")
  }
  
  if (refMeasure == 1) {
    xi <- (x-min(x))/(max(x)-min(x))
    yi <- (y-min(y))/(max(y)-min(y))
  } else {
    xi <- (x-mean(x))/(sd(x))
    yi <- (y-mean(y))/(sd(y))
  }
  x.y <- data.frame(data.x = xi, data.y = yi)
  f <- 0
  if (estimator == 1) {
    x.sort <- sort(x.y$data.x)
    y.sort <- sort(x.y$data.y)
    hx <- 0
    for(i in 1:(nx-1)) {
      delta <- x.sort[i+1] - x.sort[i]
      if (delta) {
        hx <- hx + log(abs(delta))
      }
    }
    #hx <- hx/(nx-1)+psi(nx)-psi(1)
    hx <- hx/(nx-1) + digamma(nx) - digamma(1)
    
    hy <- 0
    for(i in 1:(ny-1)) {
      delta <- y.sort[i+1]- y.sort[i]
      if (delta){
        hy <- hy + log(abs(delta))
      }
    }
    #hy <- hy/(ny-1)+psi(ny)-psi(1)
    hy <- hy/(ny-1) + digamma(ny) - digamma(1)
    
    f <- hy - hx
    return(f)
    
  } else {
    a <- 0
    b <- 0
    
    indx <- sort(x.y$data.x, index.return = T)$ix
    indy <- sort(x.y$data.y, index.return = T)$ix
    x.sortx <- x.y[indx, ]$data.x
    y.sortx <- x.y[indx, ]$data.y
    x.sorty <- x.y[indy, ]$data.x
    y.sorty <- x.y[indy, ]$data.y
    for(i in 1:(nx-1)) {
      if (x.sortx[i+1]!=x.sortx[i] & y.sortx[i+1]!=y.sortx[i]){
        a <- a + log(abs(y.sortx[i+1]-y.sortx[i])/(x.sortx[i+1]-x.sortx[i]))
      }
      if (y.sorty[i+1]!=y.sorty[i] & x.sorty[i+1]!=x.sorty[i]){
        b <- b + log(abs(x.sorty[i+1]-x.sorty[i])/(y.sorty[i+1]-y.sorty[i]))
      }
    }
    
    f <- (a-b)/nx
    return(f)
  }
  
  #return(f)
}

IGCIfor2value <- function(data, sample.size = 100, times = 100, refMeasure = "g", estimator = "e") {
  if (dim(data)[2]!=2) {
    stop("More than 2 values!\n")
  }
  if (dim(data)[1]<100) {
    stop("less than 100 obervations!\n")
  }
  F1 <- 0
  F2 <- 0
  f1 <- 0
  f2 <- 0
  n <- dim(data)[1]
  for (i in 1:times) {
    ind.sample <- sample(n, sample.size, replace = T)
    data.sample <- data[ind.sample, ]
    x <- data.sample[,1]
    y <- data.sample[,2]
    f <- igci(x, y, refMeasure = refMeasure, estimator = estimator)
    if (f<0) {
      F1 <- F1 + 1
      f1 <- f1 + f
    } else if (f>0){
      F2 <- F2 + 1
      f2 <- f2 + f
    }else {
      next
    }
  }
  f.mean <- c(f1/F1, f2/F2)
  if (F1==0){
    f.mean[1] <- NA
  }
  if (F2==0){
    f.mean[2] <- NA
  }
  if (F1 > F2) {
    p.value <- 1 - F1/times
    causal <- paste0(colnames(data)[1],"->",colnames(data)[2])
  } else if (F2 > F1) {
    p.value <- 1 - F2/times
    causal <- paste0(colnames(data)[2],"->",colnames(data)[1])
  } else {
    p.value <- NA
    causal <- "not known"
  }
  return(list(f.mean = f.mean, causal.direction = causal, p.value = p.value))
}
