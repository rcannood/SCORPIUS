#' SCORPIUS: Trajectory inference from single-cell RNA sequencing data.
#'
#' SCORPIUS orders single cells with regard to an implicit timeline,
#' such as cellular development or progression over time.
#'
#' @section Foo functions:
#' The foo functions ...
#'
#' @docType package
#' @name SCORPIUS
NULL

#' @title scRNA-seq data of dendritic cell progenitors.
#'
#' @description This dataset contains the expression values of 15752 genes for 248 dendritic cell progenitors.
#' Each cell is in one of three maturation stages: MDP, CDP or PreDC. The levels of the factor in \code{sample.info}
#' are ordered according to the maturation process.
#'
#' @format A list containing two data frames, \code{expression} (248x15752) and \code{sample.info} (248x1).
#'
#' @references Schlitzer A, Sivakamasundari V, Chen J, Sumatoh HR et al. Identification of cDC1- and cDC2-committed DC progenitors reveals early lineage priming at the common DC progenitor stage in the bone marrow. Nat Immunol 2015 Jul;16(7):718-28. PMID: 26054720
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60783}
#'
#' @docType data
"ginhoux"
