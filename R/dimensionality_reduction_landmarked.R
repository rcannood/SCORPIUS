
#' Multidimensional scaling with landmarks
#'
#' @param x a numeric matrix
#' @param dist_fun the distance function to be used; must have exactly two arguments, namely dist_fun(x, y).
#' @param ndim the maximum dimension of the space which the data are to be represented in; must be in {1, 2, \ldots, n-1}.
#' @param landmark_method Must be "naive" for now. Other landmark methods will be supported in the future.
#' @param num_landmarks the number of landmarks to be selected.
#' @param rescale A logical indicating whether or not the returned space should be rescaled and centered.
#'
#' @return A list containing the reduced space of the landmarks and the complete dataset.
#' @export
#'
#' @examples
#' data(ginhoux)
#' space <- reduce_dimensionality_landmarked(
#'   ginhoux$expression,
#'   correlation_distance,
#'   num_landmarks = 200)
#' draw_trajectory_plot(space, ginhoux$sample_info$group_name)
reduce_dimensionality_landmarked <- function(x, dist_fun, ndim = 3, landmark_method = "naive", num_landmarks = 1000, rescale = T) {
  lm_out <- landmark_selection(x, dist_fun, landmark_method, num_landmarks)
  cmd_out <- cmdscale_withlandmarks(lm_out$dist_lm, lm_out$dist_2lm, ndim = ndim, rescale = rescale)
  attr(cmd_out, "landmarks") <- lm_out$ix_lm
  cmd_out
}

landmark_selection <- function(x, dist_fun, landmark_method, num_landmarks) {
  switch(
    landmark_method,
    "naive" = {
      ix_lm <- sample.int(nrow(x), num_landmarks)
      dist_lm <- dist_fun(x[ix_lm,,drop=FALSE], x[ix_lm,,drop=FALSE])
      dist_2lm <- dist_fun(x[ix_lm,,drop=FALSE], x)
      list(ix_lm = ix_lm, dist_lm = dist_lm, dist_2lm = dist_2lm)
    },
    {
      stop("landmark_method must be ", sQuote("naive"), ".")
    }
  )
}

cmdscale_withlandmarks <- function(dist_lm, dist_2lm, ndim = 3, rescale = T) {
  d <- dist_lm
  if (anyNA(d))
    stop("NA values not allowed in 'd'")

  x <- as.matrix(d^2)
  storage.mode(x) <- "double"
  if (nrow(x) != ncol(x))
    stop("distances must be result of 'dist' or a square matrix")

  rn <- rownames(x)
  rn_all <- colnames(dist_2lm)
  n <- as.integer(nrow(x))
  N <- as.integer(ncol(dist_2lm))

  if((ndim <- as.integer(ndim)) > n - 1 || ndim < 1)
    stop("'ndim' must be in {1, 2, ..  n - 1}")

  # double center data
  mu_n <- rowMeans(x)
  mu <- mean(x)
  x_dc <- x - rep(mu_n, n) - rep(mu_n, each = n) + mu

  # classical MDS on landmarks
  e <- eigen(-x_dc/2, symmetric = TRUE)
  ev <- e$values[seq_len(ndim)]
  evec <- e$vectors[, seq_len(ndim), drop = FALSE]
  ndim1 <- sum(ev > 0)
  if (ndim1 < ndim) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", ndim1, ndim), domain = NA)
    evec <- evec[, ev > 0,  drop = FALSE]
    ev <- ev[ev > 0]
  }
  Slm <- evec * rep(sqrt(ev), each=n)

  # distance-based triangulation
  points_inv <- evec / rep(sqrt(ev), each=n)
  S <- (-t(dist_2lm - rep(mu_n, each = N))/2) %*% points_inv

  # clean up dimension names
  dimnames(Slm) <- list(rn, paste0("Comp", seq_len(ndim)))
  dimnames(S) <- list(rn_all, paste0("Comp", seq_len(ndim)))

  # rescale if necessary
  if (rescale) {
    Slm <- scale_uniform(Slm)
    S <- scale_uniform(S)
  }

  # output
  attr(S, "landmark_space") <- Slm

  S
}
