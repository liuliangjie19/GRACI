#' @title impute the SnpMatrix
#'
#' @description
#' \code{impute.res} impute the SnpMatrix and remove the remain missing value in the matrix.
#'
#' @details
#' This function mimics a Leave-One-Out process where missing SNPs are 
#' imputed for an individual based on a model trained on all other individuals.
#' This function use the output of the function \code{read.ped.to.snpMatrix}. 
#' 
#' @param res output of the function \code{read.ped.to.snpMatrix}.
#' @param removed (optional) a character string matching one of the following items: snps, sample, none. Describes the action taken in case of remaining missing values.
#' @param silence logical (default:TRUE). If \code{TRUE}, this function won't print the detail of imputation.
#' @export
#' @return a list of class GGIGCI.data with 3 items:
#' \item{SnpMatrix}{SnpMatrix contain each sample's each snp's genotype}
#' \item{gene.info}{genetic information of each snp loci}
#' \item{sample.info}{sample information contain Sex and phenotype}
#' @author Liangjie Liu <liuliangjie@@sjtu.edu.cn>
#' @examples
#' temp.im <- impute.res(res, removed = "sample")
#' 
#' 

impute.res <- function(res, removed = c("snps", "sample", "none"), silence = TRUE){
  out <- list()
  if (class(res)!="list") {
    stop("please input vaild data.\nThe data of res should be the output of the function:read.ped.to.snpMatrix().\n")
  }

  removed <- match.arg(removed)
  snp.matrix <- res$ped$genotypes
  gene.snp.info <- res$info
  #print(remain)
  #print(colnames(snp.matrix))
  #print(head(gene.snp.info))
  imputed.snp.matrix <- as(snp.matrix, "numeric")

  bar <- txtProgressBar(0, nrow(snp.matrix), char = ">", style = 3)
  for(i in 1:nrow(snp.matrix)) {
    if(all(is.na(snp.matrix[i, ]))){
      warning(paste0("the snp of ", rownames(snp.matrix)[i], "is all NA!!!"))
      next
    }
    snp.miss.idx <- which(is.na(snp.matrix[i, ]))
    if (length(snp.miss.idx)>0) {
      m <- snp.matrix[-i, snp.miss.idx]
      p <- snp.matrix[-i, -snp.miss.idx]
      #temp.snp.list <- colnames(snp.matrix)[snp.miss.idx]
      pm <- gene.snp.info[snp.miss.idx, 4]
      pp <- gene.snp.info[-snp.miss.idx, 4]
      if(silence){
        #message(as.character(i))
        r <- silence(snpStats::snp.imputation)(p, m, pp, pm)
      }
      else {
        #r <- silence(snpStats::snp.imputation)(p, m, pp, pm)
        r <- snpStats::snp.imputation(p, m, pp, pm)
      }
      temp <- snpStats::impute.snps(r, snp.matrix[i, ])

      imputed.snp.matrix[i, is.na(imputed.snp.matrix[i, ])] <- round(temp)
    }

    setTxtProgressBar(bar, i)
  }
  #out <- as(imputed.snp.matrix, "SnpMatrix")
  message("\n")

  if (any(is.na(imputed.snp.matrix))) {
    if (removed=="snps") {
      #print(remain)
      c <- which(colSums(is.na(imputed.snp.matrix))>0)
      imputed.snp.matrix <- imputed.snp.matrix[, -c]
      warning(paste0("Removed ", length(c), " SNP(s) due to remain missing value."))

      imputed.snp.matrix <- as(imputed.snp.matrix, "SnpMatrix")
      #out$ped$genotypes <- imputed.snp.matrix
      #out$info <- gene.snp.info
    }
    else if (removed=="sample") {
      #print(remain)
      r <- which(rowSums(is.na(imputed.snp.matrix))>0)
      imputed.snp.matrix <- imputed.snp.matrix[-r, ]
      warning(paste0("Removed ", length(r), " sample(s) due to remain missing value."))

      imputed.snp.matrix <- as(imputed.snp.matrix, "SnpMatrix")
      #out$ped$genotypes <- imputed.snp.matrix
    }
    else if (removed=="none") {
      #print(remain)
      warning(paste0(sum(is.na(imputed.snp.matrix)), " missing values are still in the SnpMatrix."))
    }
  }
  sample.select <- which(rownames(res$ped$fam)%in%rownames(imputed.snp.matrix))
  sample.fam.info <- res$ped$fam[sample.select, ]
  snp.select <- which(gene.snp.info[, 3]%in%colnames(imputed.snp.matrix))
  gene.snp.info <- gene.snp.info[snp.select, ]

  out <- list(SnpMatrix = imputed.snp.matrix, gene.info = gene.snp.info, sample.info = sample.fam.info)
  class(out) <- c("GGIGCI.data", "list")
  return(out)
}

silence <- function(f){
  return(function(...) {capture.output(w<-f(...));return(w);});
}
