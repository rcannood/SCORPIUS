check_numeric_matrix <- function(x, param_name, is_nullable) {
  check <- is.matrix(x) || is.data.frame(x) || (is.nullable && is.null(x))
  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be ",
      ifelse(is_nullable, "NULL, ", ""),
      "a numeric matrix, or a data frame with only numeric columns."
    )
    stop(error)
  }
}
