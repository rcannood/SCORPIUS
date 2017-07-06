
#' Multidimensional scaling with landmarks
#'
#' @param x a numeric matrix
#' @param dist.fun the distance function to be used; must have exactly two arguments, namely dist.fun(x, y).
#' @param ndim the maximum dimension of the space which the data are to be represented in; must be in {1, 2, \ldots, n-1}.
#' @param landmark.method Must be "naive" for now. Other landmark methods will be supported in the future.
#' @param num.landmarks the number of landmarks to be selected.
#' @param rescale A logical indicating whether or not the returned space should be rescaled and centered.
#'
#'
#' @return A list containing the reduced space of the landmarks and the complete dataset.
#' @export
#'
#' @examples
#' data(ginhoux)
#' space <- mds_withlandmarks(ginhoux$expression, correlation.distance)
#' draw.trajectory.plot(space, ginhoux$sample.info$group.name)
reduce.dimensionality.landmarked <- function(x, dist.fun, ndim = 3, landmark.method = "naive", num.landmarks = 1000, rescale = T) {
  lm.out <- landmark.selection(x, dist.fun, landmark.method, num.landmarks)
  cmd.out <- cmdscale.withlandmarks(lm.out$dist.lm, lm.out$dist.2lm, ndim = ndim)
  attr(cmd.out, "landmarks") <- lm.out$ix.lm
  cmd.out
}

landmark.selection <- function(x, dist.fun, landmark.method, num.landmarks) {
  switch(
    landmark.method,
    "naive" = {
      ix.lm <- sample.int(nrow(x), num.landmarks)
      dist.lm <- dist.fun(x[ix.lm,,drop=F], x[ix.lm,,drop=F])
      dist.2lm <- dist.fun(x[ix.lm,,drop=F], x)
      list(ix.lm = ix.lm, dist.lm = dist.lm, dist.2lm = dist.2lm)
    },
    {
      stop("landmark.method must be ", sQuote("naive"), "; other landmark methods will be supported in the future.")
    }
  )
}

cmdscale.withlandmarks <- function(dist.lm, dist.2lm, ndim = 3, rescale = T) {
  d <- dist.lm
  if (anyNA(d))
    stop("NA values not allowed in 'd'")

  x <- as.matrix(d^2)
  storage.mode(x) <- "double"
  if (nrow(x) != ncol(x))
    stop("distances must be result of 'dist' or a square matrix")

  rn <- rownames(x)
  rn.all <- colnames(dist.2lm)
  n <- as.integer(nrow(x))
  N <- as.integer(ncol(dist.2lm))

  if((ndim <- as.integer(ndim)) > n - 1 || ndim < 1)
    stop("'ndim' must be in {1, 2, ..  n - 1}")

  # double center data
  mu.n <- rowMeans(x)
  mu <- mean(x)
  x.dc <- x - rep(mu.n, n) - rep(mu.n, each = n) + mu

  # classical MDS on landmarks
  e <- eigen(-x.dc/2, symmetric = TRUE)
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
  points.inv <- evec / rep(sqrt(ev), each=n)
  S <- (-t(dist.2lm - rep(mu.n, each = N))/2) %*% points.inv

  # clean up dimension names
  dimnames(Slm) <- list(rn, paste0("Comp", seq_len(ndim)))
  dimnames(S) <- list(rn.all, paste0("Comp", seq_len(ndim)))

  # rescale if necessary
  if (rescale) {
    Slm <- rescale.and.center(Slm)
    S <- apply.quant.scale(S, attr(Slm, "center"), attr(Slm, "scale"))
  }

  # output
  attr(S, "landmark_space") <- Slm

  S
}
