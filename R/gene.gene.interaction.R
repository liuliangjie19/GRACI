#' @title calculate the gene-gene interaction p.val and s.val and return the result with the input data.
#'
#' @description
#' \code{gene.gene.interaction} calculate the gene-gene interaction p.val and s.val based on the method GBIGM.
#'
#' @details
#' This function use the output of the function \code{impute.res}. The 
#' class of the input data should be GGIGCI.data. The interaction calculated 
#' by this function is based on the information gained method, so the interaction 
#' can be non-linear. Larger n.times you set, more time this function run. In order 
#' to reduce the time usage, you can set the param:gene.level as minP to make this
#' function run faster. And when your input data size is too big to get a desired 
#' result, you can also try to set gene.level as minP. On the other hand, if the number
#' of SNPs on your genes are too few to present the genes, you should set gene.level as 
#' gene.level to prevent the wrong result.
#' 
#' @param data output of the function \code{impute.res}.
#' @param n.times the times of the permutation.
#' @param selected.gene.list a vector.you can selected part of genes in the gene.info.
#' @param gene.level a character string matching one of the following item: gene.level, default, minP. more information is available in the details.
#' @export
#' @return a list of class GGIUD contain 7 items:
#' \item{SnpMatrix}{function output}
#' \item{gene.info}{genetic information of each snp loci}
#' \item{sample.info}{sample information contain Sex and phenotype}
#' \item{GGI.ud}{gene-gene interaction table, you can plot a network by using this table and the package\code{igrah}}
#' \item{GENE.table}{the selected gene list and the number of snps on each gene}
#' \item{P.matrix}{p value matrix}
#' \item{S.matrix}{delta H matrix}
#' @author Liangjie Liu <liuliangjie@@sjtu.edu.cn>
#' @examples
#' temp.ggiud <- gene.gene.interaction(temp.im, n.times = 1000, selected.gene.list = c("DRD1", "DRD2", "DRD3"))
#'
#'

gene.gene.interaction <- function(data, n.times = 1000, selected.gene.list=NULL, gene.level = c("gene.level", "default", "minP")){

  if (!inherits(data, "GGIGCI.data")) {
    stop("data must be one of the class 'GGIGCI.data'(the output of the function impute.res).")
  }
  snp.info <- data$SnpMatrix
  gene.info <- data$gene.info
  sample.info <- data$sample.info
  gene.level <- match.arg(gene.level)
  #phen <- sample.info[,6]
  if (is.null(selected.gene.list)) {
    selected.gene.list <- unique(gene.info[, 2])
  }
  else if (any(!selected.gene.list%in%gene.info[,2])) {
    stop("invaild gene selected.")
  }
  #else {
  gene.info <- gene.info[which(gene.info[, 2]%in%selected.gene.list), ]
  #}
  if (any(sample.info[,2]!=rownames(snp.info))){
    warning("miss match between sample info matrix and SnpMatrix.")
  }
  phen <- sample.info[,6]
  p.val.matrix <- diag(0, nrow = length(selected.gene.list))
  colnames(p.val.matrix) <- selected.gene.list
  rownames(p.val.matrix) <- selected.gene.list
  s.val.matrix <- diag(0, nrow = length(selected.gene.list))
  colnames(s.val.matrix) <- selected.gene.list
  rownames(s.val.matrix) <- selected.gene.list

  interaction.pairs <- combn(selected.gene.list, 2)
  gene.tabe <- data.frame(Genenames = selected.gene.list, Snps = table(gene.info[,2])[selected.gene.list])[,-2]
  ggi.ud <- data.frame(from = interaction.pairs[1, ], to = interaction.pairs[2, ],
                       P.val = rep(0, length(interaction.pairs[1, ])), S.val = rep(0, length(interaction.pairs[1, ])))
  nc.i <- ncol(interaction.pairs)

  bar <- txtProgressBar(0, nc.i, char = ">", style = 3)
  for(n in 1:nc.i){
    gene1 <- interaction.pairs[1, n]
    gene2 <- interaction.pairs[2, n]
    #print(gene1)
    #print(gene2)
    snp.gene1 <- gene.info[which(gene.info[, 2]==gene1), 3]
    snp.gene2 <- gene.info[which(gene.info[, 2]==gene2), 3]
    snp.info.gene1 <- snp.info[, colnames(snp.info)%in%snp.gene1]
    snp.info.gene2 <- snp.info[, colnames(snp.info)%in%snp.gene2]
    #temp.result <- ifelse(gene.level, GBIGM.gene.level(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times),
                          #GBIGM(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times))
    if (gene.level == "gene.level") {
      temp.result <- GBIGM.gene.level(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times)
    }
    else if (gene.level == "minP") {
      #snp.info.gene1 <- fex(snp.info.gene1, gene.info[which(gene.info[, 2]==gene1), ])
      #snp.info.gene2 <- fex(snp.info.gene2, gene.info[which(gene.info[, 2]==gene2), ])
      division.gene1 <- cluster.gene(snp.info.gene1)
      division.gene2 <- cluster.gene(snp.info.gene2)

      pair.start <- expand.grid(division.gene1$start, division.gene2$start)
      pair.end <- expand.grid(division.gene1$end, division.gene2$end)

      temp.pair <- data.frame(start.gene1 = pair.start[,1], end.gene1 = pair.end[, 1], start.gene2 = pair.start[,2], end.gene2 = pair.end[, 2])
      p.list <- rep(NA, nrow(temp.pair))
      s.list <- rep(NA, nrow(temp.pair))

      for(i in 1:nrow(temp.pair)){
        sub.gene1 <- (temp.pair$start.gene1[i]):(temp.pair$end.gene1[i])
        sub.gene2 <- (temp.pair$start.gene2[i]):(temp.pair$end.gene2[i])
        sub.p.s <- GBIGM(phen = phen, snp.info.gene1[, sub.gene1], snp.info.gene2[, sub.gene2], n.times = n.times)
        p.list[i] <- sub.p.s[1]
        s.list[i] <- sub.p.s[2]
      }
      
      p.list <- p.adjust(p.list, "BH")
      temp.result <- c(min(p.list), s.list[which.min(p.list)])
    }
    else if (gene.level == "default") {
      temp.result <- GBIGM(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times)
      #temp.result <- GBIGM.test(Y = phen, G1 = snp.info.gene1, G2 = snp.info.gene2, n.perm = n.times)
    }
    p.val.matrix[gene1, gene2] <- p.val.matrix[gene2, gene1] <- temp.result[1]
    #p.val.matrix[gene1, gene2] <- p.val.matrix[gene2, gene1] <- temp.result$p.val
    s.val.matrix[gene1, gene2] <- s.val.matrix[gene2, gene1] <- temp.result[2]
    #s.val.matrix[gene1, gene2] <- s.val.matrix[gene2, gene1] <- temp.result$statistic
    #ggi.ud[n, 3:4] <- c(temp.result$p.val, temp.result$statistic)
    ggi.ud[n, 3:4] <- temp.result

    setTxtProgressBar(bar, n)
  }

  message("\n")
  snp.matrix <- snp.info[, colnames(snp.info)%in%gene.info[,3]]
  out <- list(SnpMatrix = snp.matrix, gene.info = gene.info, sample.info = sample.info,
              GGI.ud = ggi.ud, GENE.table = gene.tabe, P.matrix = p.val.matrix, S.matrix = s.val.matrix)
  class(out) <- c("GGIUD", "list")
  return(out)
}

cluster.gene <- function(gene){
  if (ncol(gene)>20){
    distance.gene <- snpStats::ld(gene, gene, stats = "R.squared")
    distance.gene <- as.dist(1 - distance.gene)
    clust.tree.gene <- rioja::chclust(distance.gene)
    k.gene <- cutree(clust.tree.gene, k = 1:(ncol(gene)-19))
    max.gene <- sapply(1:(ncol(gene)-19),FUN=function(i){return(max(table(as.factor(k.gene[,i]))))})
    id.gene <- which(max.gene<=20)[1]
  }
  else if (ncol(gene)==1){
    return(list(start = c(1), end = c(1)))
  }
  else {
    k.gene <- matrix(rep(1, ncol(gene)), ncol = 1)
    rownames(k.gene) <- colnames(gene)
    max.gene <- ncol(gene)
    id.gene <- 1
  }
  #print(sapply(1:(length(k.gene[,id.gene])-1),FUN=function(i){return(k.gene[,id.gene][i]!=k.gene[,id.gene][i+1])}))

  division.gene.start <- c(1,1+as.numeric(which(sapply(1:(length(k.gene[,id.gene])-1),FUN=function(i){return(k.gene[,id.gene][i]!=k.gene[,id.gene][i+1])}))))
  division.gene.end <- c(as.numeric(which(sapply(1:(length(k.gene[,id.gene])-1),FUN=function(i){return(k.gene[,id.gene][i]!=k.gene[,id.gene][i+1])}))),ncol(gene))

  return(list(start = division.gene.start, end = division.gene.end))
}

GBIGM <- function(phen, gene1, gene2, n.times = 1000) {
  x1 <- as(gene1, "numeric")
  x2 <- as(gene2, "numeric")
  #selfhead(x1)
  #selfhead(x2)
  if (nrow(x1)!=nrow(x2)){
    stop("The two SnpMatrix must have the same row numbers. ")
  }
  else if (length(phen)!=nrow(x1)){
    stop("The length of phen must be equal to the SnpMatrix's row.number.")
  }
  num.sample <- length(phen)
  
  g1 <- apply(x1, 1, paste, collapse = " ")
  g2 <- apply(x2, 1, paste, collapse = " ")
  g1.2 <- apply(cbind(x1, x2), 1, paste, collapse = " ")
  
  gx1 <- apply(cbind(phen, x1), 1, paste, collapse = " ")
  gx2 <- apply(cbind(phen, x2), 1, paste, collapse = " ")
  gx1.2 <- apply(cbind(phen, cbind(x1, x2)), 1, paste, collapse = " ")
  #print(table(gx1))
  #print(table(gx2))
  #print(table(gx1.2))
  HG1 <- entropy.vec(g1)
  HG2 <- entropy.vec(g2)
  HG1.2 <- entropy.vec(g1.2)
  
  HGx1 <- entropy.vec(gx1)-HG1
  HGx2 <- entropy.vec(gx2)-HG2
  HGx1.2 <- entropy.vec(gx1.2)-HG1.2
  
  if (min(HGx1, HGx2)==0){
    Delta1.2init <- min(HGx1, HGx2)-HGx1.2
  }
  else {
    Delta1.2init <- (min(HGx1, HGx2)-HGx1.2)/min(HGx1, HGx2)
  }
  
  D <- list()
  for(i in seq_len(n.times)){
    D[[i]] <- phen[sample(1:num.sample, num.sample)]
  }
  Delta1.2 <- rep(NA, n.times)
  for (i in seq_len(n.times)){
    #D[[i]] <- sample(1:num.sample, num.sample)
    #gs1 <- gx1[D[[i]]]
    gs1 <- apply(cbind(D[[i]], x1), 1, paste, collapse = " ")
    #gs2 <- gx2[D[[i]]]
    gs2 <- apply(cbind(D[[i]], x2), 1, paste, collapse = " ")
    #gs1.2 <- gx1.2[D[[i]]]
    gs1.2 <- apply(cbind(D[[i]], cbind(x1, x2)), 1, paste, collapse = " ")
    
    HGs1 <- entropy.vec(gs1)-HG1
    HGs2 <- entropy.vec(gs2)-HG2
    HGs1.2 <- entropy.vec(gs1.2)-HG1.2
    if (min(HGs1, HGs2)==0){
      Delta1.2tmp <- min(HGs1, HGs2)-HGs1.2
    }
    else {
      Delta1.2tmp <- (min(HGs1, HGs2)-HGs1.2)/min(HGs1, HGs2)
    }
    Delta1.2[i] <- Delta1.2tmp
  }
  pval <- sum(Delta1.2>=Delta1.2init)/n.times
  
  stat <- Delta1.2init
  
  return(c(pval, stat))
}

###when the number of snp loci on the your gene is too small, the gene.level.method can be considered.
GBIGM.gene.level <- function(phen, gene1, gene2, n.times = 1000) {
  
  x1 <- as(gene1, "numeric")
  x2 <- as(gene2, "numeric")
  xt1 <- apply(x1, 2, correct.genotype)
  xt2 <- apply(x2, 2, correct.genotype)
  
  if (nrow(x1)!=nrow(x2)){
    stop("Error 1.")
  }
  else if (length(phen)!=nrow(x1)){
    stop("Error 2")
  }
  num.sample <- length(phen)
  
  g1 <- apply(xt1, 1, is.mutation)
  g2 <- apply(xt2, 1, is.mutation)
  #g1.2 <- apply(cbind(xt1, xt2), 1, is.mutation)
  g1.2 <- apply(cbind(g1, g2), 1, paste, collapse = " ")
  
  gx1 <- apply(cbind(phen, g1), 1, paste, collapse = " ")
  gx2 <- apply(cbind(phen, g2), 1, paste, collapse = " ")
  gx1.2 <- apply(cbind(phen, g1.2), 1, paste, collapse = " ")
  
  HG1 <- entropy.vec(g1)
  HG2 <- entropy.vec(g2)
  HG1.2 <- entropy.vec(g1.2)
  
  HGx1 <- entropy.vec(gx1)-HG1
  HGx2 <- entropy.vec(gx2)-HG2
  HGx1.2 <- entropy.vec(gx1.2)-HG1.2
  
  if (min(HGx1, HGx2)==0){
    Delta1.2init <- min(HGx1, HGx2)-HGx1.2
  }
  else {
    Delta1.2init <- (min(HGx1, HGx2)-HGx1.2)/min(HGx1, HGx2)
  }
  
  D <- list()
  for(i in seq_len(n.times)){
    D[[i]] <- phen[sample(1:num.sample, num.sample)]
  }
  Delta1.2 <- rep(NA, n.times)
  for (i in seq_len(n.times)){
    #D[[i]] <- sample(1:num.sample, num.sample)
    #gs1 <- gx1[D[[i]]]
    gs1 <- apply(cbind(D[[i]], g1), 1, paste, collapse = " ")
    #gs2 <- gx2[D[[i]]]
    gs2 <- apply(cbind(D[[1]], g2), 1, paste, collapse = " ")
    #gs1.2 <- gx1.2[D[[i]]]
    gs1.2 <- apply(cbind(D[[i]], g1.2), 1, paste, collapse = " ")
    
    HGs1 <- entropy.vec(gs1)-HG1
    HGs2 <- entropy.vec(gs2)-HG2
    HGs1.2 <- entropy.vec(gs1.2)-HG1.2
    if (min(HGs1, HGs2)==0){
      Delta1.2tmp <- min(HGs1, HGs2)-HGs1.2
    }
    else {
      Delta1.2tmp <- (min(HGs1, HGs2)-HGs1.2)/min(HGs1, HGs2)
    }
    Delta1.2[i] <- Delta1.2tmp
  }
  pval <- sum(Delta1.2>=Delta1.2init)/n.times
  
  stat <- Delta1.2init
  
  return(c(pval, stat))
}

correct.genotype <- function(snp.x){
  w <- "Wildtype"
  m <- "Mutation"
  
  #snp.table <- table(snp.x)
  l0 <- length(snp.x[which(snp.x==0)])
  l2 <- length(snp.x[which(snp.x==2)])
  #if (snp.table[[1]]<snp.table[[2]]){
  if (l0<l2){
    snp.t <- ifelse(snp.x==2, w, m)
  }
  else {
    snp.t <- ifelse(snp.x==0, w, m)
  }
  
  return(snp.t)
}

is.mutation <- function(xt){
  return(ifelse("Wildtype"%in%xt, 1, 0))
}

entropy.vec <- function(X){
  pourc<-percentofX(X)
  H<-(-sum(pourc*log(pourc)))
  return(H)
}

percentofX <- function(X){
  n.ind<-length(X)
  b10<-as.factor(X)
  vect<-tabulate(b10)
  pourc<-vect/n.ind
  return(pourc)
}
