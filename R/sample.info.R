#' @title sample.info data for 1536 samples
#'
#' @description A dataset containing the information of 1536 samples and their basic information.
#'
#' @docType data
#' @keywords datasets
#' @name sample.info
#' @usage sample.info
#' @format A data frame with 1536 rows and 6 variables. The variables
#' are as follows:
#' \describe{
#'   \item{Fam.id}{the family ID of each sample, if not known, family id is same to the Sample.Name}
#'   \item{Sample.Name}{The unique ID of each sample}
#'   \item{Paternal.ID}{The paternal ID of each sample, if not known, paternal id = 0}
#'   \item{Maternal.ID}{The maternal ID of each sample, if not known, maternal id = 0}
#'   \item{Sex}{The sex of each sample, 1 for male; 2 for female; other for unknown}
#'   \item{Phen}{The Phenotype of each sample, 0=unknown; 1=unaffected; 2=affected}
#' }
#' @source unknown.
NULL