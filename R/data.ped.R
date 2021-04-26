#' @title ped data for 1536 samples and 92 snps
#'
#' @description A dataset containing the information of 1536 samples
#' and their 92 snp locis data.
#'
#' @docType data
#' @keywords datasets
#' @name data.ped
#' @usage data.ped
#' @format A data frame with 1536 rows and 98 variables. The variables
#' are as follows:
#' \describe{
#'   \item{Well}{The Well number in each plate}
#'   \item{Sample.Name}{The ID of each sample}
#'   \item{Paternal.ID}{The ID of the father of each sample, 0 for unkonw}
#'   \item{Maternal.ID}{The ID of the mother of each sample, 0 for unkonw}
#'   \item{Sex}{The sex status of each sample;1 for male; 2 for female; other for unknown}
#'   \item{Phen}{The Phenotype of each sample, 0=unknown; 1=unaffected; 2=affected}
#'   \item{snp locis}{The ID of each snp loci in the Well plate}
#' }
#' @source unknown.
NULL
