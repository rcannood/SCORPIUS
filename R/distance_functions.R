#' Deprecated distance functions
#'
#' Passing `euclidean_distance()` and `correlation_distance()` to [reduce_dimensionality()] is deprecated.
#' Instead simply pass `"euclidean"` or `"pearson"`, respectively.
#'
#' @param x deprecated
#' @param y deprecated
#'
#' @rdname distance
#' @export
euclidean_distance <- function(x, y) {
  .Deprecated(new = "\"euclidean\"")
}

#' @rdname distance
#' @export
correlation_distance <- function(x, y) {
  .Deprecated(new = "\"pearson\"")
}
