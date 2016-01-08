# #' Generate R documentation from inline comments.
# #'
# #' Roxygen2 allows you to write documentation in comment blocks co-located
# #' with code.
# #'
# #' The only function you're likely to need from \pkg{roxygen2} is
# #' \code{\link{roxygenize}}. Otherwise refer to the vignettes to see
# #' how to format the documentation.
# "_PACKAGE"


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
"ginhoux"