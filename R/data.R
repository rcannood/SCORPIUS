#' @title Generate a synthetic dataset
#'
#' @description \code{generate_dataset} generates an synthetic dataset which can be used for visualisation purposes.
#'
#' @usage
#' generate_dataset(
#'   num_samples = 400,
#'   num_genes = 500,
#'   num_groups = 4
#' )
#'
#' @param num_samples The number of samples the dataset will contain.
#' @param num_genes The number of genes the dataset will contain.
#' @param num_groups The number of groups the samples will be split up in.
#'
#' @return A list containing the expression data and the meta data of the samples.
#'
#' @importFrom stats poly rnorm runif
#'
#' @export
#'
#' @seealso [SCORPIUS]
#'
#' @examples
#' ## Generate a dataset
#' dataset <- generate_dataset(num_genes = 500, num_samples = 1000, num_groups = 4)
#'
#' ## Reduce dimensionality and infer trajectory with SCORPIUS
#' space <- reduce_dimensionality(dataset$expression, ndim = 2)
#' traj <- infer_trajectory(space)
#'
#' ## Visualise
#' draw_trajectory_plot(space, path=traj$path, progression_group=dataset$sample_info$group_name)
generate_dataset <- function(num_samples = 400, num_genes = 500, num_groups = 4) {
  # make names for each group, gene and sample
  group_names <- paste0("Group ", seq_len(num_groups))
  gene_names <- paste0("Gene", seq_len(num_genes))
  sample_names <- paste0("Sample", seq_len(num_samples))


  # construct the sample info
  x <- seq(-1, 1, length.out = num_samples)
  group <- cut(x, breaks = num_groups, labels = group_names)
  sample_info <- data.frame(row.names = sample_names, group_name = group)

  # apply function and determine noise sd
  y <- stats::poly(x, 2)
  sd <- .012 * sqrt(num_genes)

  # generate expression data
  expression <- sapply(seq_len(num_genes), function(g) {
    scale <- stats::rnorm(ncol(y), mean=0, sd=1)
    noise <- stats::rnorm(length(x), sd=sd)
    rowSums(sweep(y, 2, scale, "*")) + noise
  })
  dimnames(expression) <- list(sample_names, gene_names)

  # simulate genes that are not expressed
  weighted_random_sample <- function(data, weights, n){
    key <- stats::runif(length(data)) ^ (1 / weights)
    data[order(key, decreasing=TRUE)][seq_len(n)]
  }

  undetectable <- which(expression < 0)
  undetectable <- weighted_random_sample(undetectable, -expression[undetectable], round(length(undetectable)*.5))

  # shift expression
  expression <- expression + .5

  # set everything below 0 or that was marked as undetectable to zero
  expression[expression < 0 | seq_along(expression) %in% undetectable] <- 0

  # rescale to a reasonable value
  expression <- expression / max(expression) * 20

  list(expression = expression, sample_info = sample_info)
}

#' @title scRNA-seq data of dendritic cell progenitors.
#'
#' @description This dataset contains the expression values of the
#' top 2000 most variable genes for 248 dendritic cell progenitors.
#' Each cell is in one of three maturation stages: MDP, CDP or PreDC.
#' The levels of the factor in \code{sample.info} are ordered according to the maturation process.
#'
#' The number of genes had to be reduced specifically for reducing the package size of SCORPIUS.
#' Use the following code to download the original data:
#' \preformatted{
#' download.file("https://github.com/rcannood/SCORPIUS/raw/master/data-raw/ginhoux_orig.rds", destfile = "local.rds")
#' ginhoux <- readRDS("local.rds")
#' # do something with ginhoux
#' }
#'
#' @format A list containing two data frames, \code{expression} (248x2000) and \code{sample_info} (248x1).
#'
#' @references Schlitzer A, Sivakamasundari V, Chen J, Sumatoh HR et al.
#' Identification of cDC1- and cDC2-committed DC progenitors reveals early lineage priming at
#'  the common DC progenitor stage in the bone marrow. Nat Immunol 2015 Jul;16(7):718-28. PMID: 26054720
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60783}
#'
#' @docType data
#'
#' @seealso [SCORPIUS]
"ginhoux"
