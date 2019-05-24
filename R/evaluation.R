#' Deprecated evaluation functions
#'
#' These functions are outdated. Use \code{dyneval::evaluate_ti_method} instead. For more info, visit \url{https://dynverse.org}.
#'
#' @param time deprecated
#' @param progression deprecated
#'
#' @rdname deprecated
#'
#' @export
evaluate_trajectory <- function(time, progression) {
  .Deprecated(new = "evaluate_ti_method", package = "dyneval")
}

#' @rdname deprecated
#'
#' @param space deprecated
#' @param k deprecated
#'
#' @export
evaluate_dim_red <- function(space, progression, k = 5) {
  .Deprecated(new = "evaluate_ti_method", package = "dyneval")
}
