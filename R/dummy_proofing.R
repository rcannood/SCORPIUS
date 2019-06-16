#' @importFrom dynutils is_sparse
check_numeric_matrix <- function(x, param_name, is_nullable = FALSE, finite = FALSE, sparse = FALSE) {
  if (is_nullable && is.null(x)) {
    return(invisible())
  }

  check <- is.matrix(x) || is.data.frame(x)
  if (sparse) {
    check <- check || dynutils::is_sparse(x)
  }

  if (check) {
    for (j in seq_len(ncol(x))) {
      check <- check && is.numeric(x[,j]) && (!finite || all(is.finite(x[,j])))
    }
  }

  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be ",
      ifelse(is_nullable, "NULL, ", ""),
      "a numeric matrix, ",
      ifelse(sparse, "a sparse numeric matrix, ", ""),
      " or a data frame containing only ",
      ifelse(finite, "finite ", ""),
      "numeric values."
    )
    stop(error)
  }
}

check_numeric_vector <- function(x, param_name, is_nullable = TRUE, finite = FALSE, whole = FALSE, range = NULL, length = NULL, factor = FALSE) {
  if (is_nullable && is.null(x)) {
    return(invisible())
  }

  if (factor && is.factor(x)) {
    x <- as.numeric(x)
  }

  check <- is.numeric(x)

  check <- check && (!finite || all(is.finite(x)))
  check <- check && (!whole || all(round(x) == x))
  check <- check && (is.null(range) || all(range[[1]] <= x & x <= range[[2]]))
  check <- check && (is.null(length) || length(x) == length)

  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be a numeric vector consisting of ",
      ifelse(!is.null(length), paste0(length, " "), ""),
      ifelse(finite, "finite ", ""),
      ifelse(whole, "whole ", ""),
      "number(s)",
      ifelse(!is.null(range), paste0(" within the range of [", range[[1]], ", ", range[[2]], "]"), "")
    )

    stop(error)
  }
}


check_logical_vector <- function(x, param_name, is_nullable = TRUE, length = NULL) {
  if (is_nullable && is.null(x)) {
    return(invisible())
  }

  check <- is.logical(x)

  check <- check && (is.null(length) || length(x) == length)

  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be a logical vector consisting of ",
      ifelse(!is.null(length), paste0(length, " "), ""),
      "logical(s)"
    )

    stop(error)
  }
}
