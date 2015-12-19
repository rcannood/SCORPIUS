euclidean.distance <- function (x, y) {
  if (ncol(x) != ncol(y)) 
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (k in 1:nrow(y)) {
    z[, k] <- sqrt(colSums((t(x) - y[k, ])^2))
  }
  z
}

knn.distances <- function(dists, k) {
  t(apply(dists, 1, function(x) sort(x)[seq_len(min(k, length(x)))]))
}

#' Title
#'
#' @param counts 
#'
#' @return
#' @export
#'
#' @examples
calculate.distance <- function(counts) {
  cor <- cor(t(counts), method="spearman")
  sim <- (cor+1)/2
  dist <- 1 - sim
  dist
}

#' Title
#'
#' @param dist 
#' @param max.rem.pct 
#'
#' @return
#' @export
#'
#' @importFrom fitdistrplus fitdist
#' @examples
calculate.outlier <- function(dist, max.rem.pct=.9) {
  requireNamespace("fitdistrplus")
  rem <- numeric(0)
  results <- list()
  knn10 <- function(x) { sort(x, decreasing = F)[seq_len(min(10, length(x)))] }
  while (length(rem) <= max.rem.pct*(nrow(dist))+1) {
    filt <- seq_len(nrow(dist)) %in% rem
    dist.agg <- sapply(which(!filt), function(i) {
      mean(knn10(dist[i,!filt & seq_len(ncol(dist)) != i]))
    })
    fit <- fitdistrplus::fitdist(dist.agg, distr="norm")
    removed <- which(!filt)[which.max(dist.agg)]
    rem <- c(rem, removed)
    results <- c(results, list(list(filt=filt, fit=fit)))
  }
  logliks <- sapply(results, function(r) r$fit$loglik)
  result <- results[[which.max(logliks)]]
  list(is.outlier=result$filt, loglikelihoods=logliks)
}

#' Title
#'
#' @param space 
#' @param center 
#' @param max.range 
#'
#' @return
#' @export
#'
#' @examples
rescale.and.center <- function(space, center=0, max.range=1) {
  mins <- apply(space, 2, min)
  maxs <- apply(space, 2, max)
  old.center <- (maxs + mins) / 2
  scale <- max(maxs - mins)
  t(apply(space, 1, function(x) (x-old.center+center) / scale))
}

#' Title
#'
#' @param dis 
#' @param ndim 
#' @param scaled 
#'
#' @return
#' @export
#'
#' @examples
reduce.dimensionality <- function(dist, ndim, scaled=T) {
  space <- cmdscale(dist, k = ndim)
  if (scaled) {
    space <- rescale.and.center(space)
  }
  colnames(space) <- paste("Comp", seq_len(ncol(space)), sep="")
  space
}

#' Title
#'
#' @param space 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
cluster <- function(space, k) {
  fit <- kmeans(space, centers = k)
  list(
    labels=fit$cluster, 
    centers=fit$centers
  )
}

#' Title
#'
#' @param centers 
#'
#' @return
#' @export
#'
#' @examples
cluster.distance <- function(centers) {
  eucl.dist <- as.matrix(dist(centers))
  k <- nrow(centers)
  density.dist <- sapply(seq_len(k), function(i) {
    sapply(seq_len(k), function(j) {
      if (i == j) {
        0
      } else {
        twocent <- centers[c(i,j),]
        unit <- twocent[2,] - twocent[1,]
        unit <- unit / sqrt(sum(unit^2)) * .01
        num.pts <- ceiling(mean((twocent[2,] - twocent[1,])/unit))
        segment.pts <- apply(twocent, 2, function(x) seq(x[[1]], x[[2]], length.out = num.pts))
        dists <- euclidean.distance(segment.pts, space)
        mean(knn.distances(dists, 10))
      }
    })
  })
  eucl.dist * density.dist
}


#' Title
#'
#' @param centers 
#'
#' @return
#' @export
#' 
#' @importFrom GA ga
#' 
#' @examples
shortest.path <- function(centers) {
  requireNamespace("GA")
  cluster.distances <- cluster.distance(centers)
  fitness <- function(ord) {
    -sum(mapply(ord[-1], ord[-length(ord)], FUN=function(i, j) cluster.distances[[i, j]]))
  }
  ga.fit <- GA::ga(type="permutation", fitness=fitness, min=1, max=nrow(centers), monitor=function(...) { })
  centers[ga.fit@solution[1,],,drop=F]
}

#' Title
#'
#' @param space 
#' @param k 
#' @param ... 
#'
#' @return
#' @export
#' 
#' @importFrom princurve principal.curve
#' @importFrom dplyr percent_rank
#'
#' @examples
infer.trajectory <- function(space, k, ...) {
  requireNamespace("princurve")
  requireNamespace("dplyr")
  clust <- cluster(space, k)
  start <- shortest.path(clust$centers)
  fit <- princurve::principal.curve(space, start=start, plot.true=F, trace=F, stretch = 0, ...)
  
  path <- fit$s[fit$tag,,drop=F]
  colnames(path) <- paste0("Comp", seq_len(ncol(path)))
  rownames(path) <- NULL
  
  value <- dplyr::percent_rank(order(fit$tag))
  names(value) <- rownames(space)
  
  list(clustering=clust, initial.path=start, final.path=path, time=value)
}

#' Title
#'
#' @param trajectory 
#'
#' @return
#' @export
#'
#' @examples
reverse.trajectory <- function(trajectory) {
  trajectory$time <- 1-trajectory$time
  trajectory$final.path <- trajectory$final.path[rev(seq_len(nrow(trajectory$final.path))),,drop=F]
  trajectory$initial.path <- trajectory$initial.path[rev(seq_len(nrow(trajectory$initial.path))),,drop=F]
  trajectory
}

#' Title
#'
#' @param counts 
#' @param time 
#' @param degrees.of.freedom 
#' @param p.adjust.method
#'
#' @return
#' @export
#' 
#' @importFrom splines ns
#' @importFrom pbapply pblapply
#' @importFrom parallel mclapply
#' @importFrom dplyr arrange
#'
#' @examples
find.trajectory.aligned.genes <- function(counts, time, degrees.of.freedom=8, p.adjust.method="BH", q.value.cutoff=1e-10, mc.cores=1) {
  lapply.fun <- 
    if (mc.cores==1) {
      function(...) pbapply::pblapply(...)
    } else {
      function(...) parallel::mclapply(..., mc.cores=mc.cores)
    }
  runs <- lapply.fun(seq_len(ncol(counts)), function(i) {
    fit <- glm(counts[,i]~splines::ns(time, degrees.of.freedom), family=gaussian(), epsilon=1e-5)
    fit1 <- glm(counts[,i]~1, family=gaussian(), epsilon=1e-5)
    test <- anova(fit1, fit, test = "F")
    pval <- test[6][2, 1]
    
    smooth <- predict(fit, splines::ns(time, degrees.of.freedom))
    list(pval=pval, smooth=smooth)
  })
  pvals <- sapply(runs, function(x) x$pval)
  pvals[!is.finite(pvals)] <- 1
  qvals <- p.adjust(pvals, method=p.adjust.method, n=ncol(counts))
  smooth.expression <- sapply(runs, function(x) x$smooth)
  dimnames(smooth.expression) <- dimnames(counts)
  
  breaks <- c(1, .05, .001, 1e-5, 1e-10, 1e-20, 1e-40, -Inf)
  qval.cat <- cut(qvals, breaks=breaks)
  levels(qval.cat) <- c("p <= 1e-40", "1e-40 < p <= 1e-20", "1e-20 < p <= 1e-10", "1e-10 < p <= 1e-5", "1e-5 < p <= 0.001", "0.001 < p <= 0.05", "p > 0.05")
  
  is.tag <- qvals < q.value.cutoff
  
  p.values <- data.frame(gene=colnames(counts), p.value=pvals, q.value=qvals, is.tag=is.tag, category=qval.cat, stringsAsFactors = F)
  p.values <- dplyr::arrange(p.values, q.value)
  
  genes <- p.values$gene[p.values$is.tag]
  list(genes=genes, p.values=p.values, smooth.expression=smooth.expression)
}

#' Title
#'
#' @param counts 
#' @param tag.genes 
#' @param time 
#'
#' @return
#' @export
#' 
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom WGCNA mergeCloseModules
#' @importFrom dplyr bind_rows arrange
#'
#' @examples
find.modules <- function(smooth.expression, tag.genes) {
  smooth.expr <- smooth.expression[,tag.genes,drop=F]
  dissim <- calculate.distance(t(smooth.expr))
  hcl <- hclust(as.dist(dissim), method="average")
  
  labels <- dynamicTreeCut::cutreeDynamic(hcl, distM=dissim, cutHeight = 0.8, deepSplit=1, pamRespectsDendro = F, method="hybrid")
  labels2 <- WGCNA::mergeCloseModules(smooth.expr, labels, cutHeight = 0.3)$colors
  labels2 <- match(labels2, unique(labels2))
  
  modules <- dplyr::bind_rows(lapply(unique(labels2), function(l) {
    ix <- which(labels2==l)
    dimred <- reduce.dimensionality(dissim[ix, ix, drop=F], ndim=2)
    value <- dimred[,1]
    data.frame(gene=tag.genes[ix], index=ix, module=l, value=value, stringsAsFactors = F)
  }))
  modules <- dplyr::arrange(modules, module, value)
  
  modules
}


#' Title
#'
#' @param space 
#' @param colour 
#' @param contour
#'
#' @return
#' @export
#' 
#' @import ggplot2
#' @importFrom MASS kde2d
#' @importFrom dplyr bind_rows
#'
#' @examples
plot.dimensionality.reduction <- function(space, colour=NULL, contour=F) {
  min <- min(space[,1:2])
  max <- max(space[,1:2])
  diff <- (max - min)/2
  
  g <- ggplot2::ggplot() + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x="Component 1", y="Component 2", colour="Group", fill="Group") + 
    ggplot2::xlim(min-diff, max+diff) + 
    ggplot2::ylim(min-diff, max+diff) + 
    ggplot2::coord_equal(xlim=c(min-0.1*diff, max+0.1*diff), ylim=c(min-0.1*diff, max+0.1*diff))

  if (contour) {
    if (is.null(colour)) {
      g <- g + ggplot2::stat_density2d(ggplot2::aes(Comp1, Comp2), breaks=1, geom="polygon", data.frame(space), alpha=.1)
    } else {
      g <- g + ggplot2::stat_density2d(ggplot2::aes(Comp1, Comp2, fill=colour), breaks=1, geom="polygon", data.frame(space, colour), alpha=.1)
    }
  }
  
  if (is.null(colour)) {
    g <- g + ggplot2::geom_point(ggplot2::aes(Comp1, Comp2), data.frame(space))
  } else {
    g <- g + ggplot2::geom_point(ggplot2::aes(Comp1, Comp2, colour=colour), data.frame(space, colour))
  }
  
  g
}

#' Title
#'
#' @param space 
#' @param path 
#' @param colour 
#'
#' @return
#' @export
#' 
#' @import ggplot2
#'
#' @examples
plot.trajectory <- function(space, path, colour=NULL, contour=F) {
  plot.dimensionality.reduction(space, colour=colour, contour=contour) + ggplot2::geom_path(ggplot2::aes(Comp1, Comp2), data.frame(path))
}

#' Title
#'
#' @param time 
#' @param progression 
#' @param colour 
#'
#' @return
#' @export
#' 
#' @import ggplot2
#'
#' @examples
plot.trajectory.density <- function(time, progression, colour=NULL) {
  ggplot2::ggplot() + ggplot2::geom_density(ggplot2::aes(time, colour=colour), data.frame(time, colour=factor(progression))) + ggplot2::theme_classic()
}

#' Title
#'
#' @param tags 
#' @param progression.str 
#' @param time 
#' @param modules 
#'
#' @return
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#'
#' @examples
plot.modules.heatmap <- function(smooth.expression, tag.genes, progression.str, time, modules) {
  smooth.expr <- smooth.expression[,tag.genes,drop=F]
  smooth.expr.ord <- smooth.expr[order(time), modules$index]
  smooth.expr.ord.scaled <- t(scale(smooth.expr.ord))
  gaps <- which(modules$module[-1] != modules$module[-length(modules$module)])
  col.ann <- data.frame(row.names = rownames(smooth.expression), Progression=progression.str, Time=time)
  row.ann <- data.frame(row.names = modules$gene, Module=paste0("Module ", modules$module))
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  num.modules <- length(unique(modules$module))
  Module.colours <- setNames(RColorBrewer::brewer.pal(max(num.modules, 3), "Set2"), paste0("Module ", seq_len(num.modules)))
  Prog.colours <- setNames(gg_color_hue(length(levels(progression.str))), levels(progression.str))
  ann.col <- list(
    Module=Module.colours,
    Time=Prog.colours,
    Progression=Prog.colours
  )
  
  pheatmap::pheatmap(
    smooth.expr.ord.scaled, 
    cluster_rows=F, 
    cluster_cols=F, 
    gaps_row=gaps, 
    annotation_row = row.ann, 
    annotation_col = col.ann, 
    annotation_colors = ann.col, 
    breaks=c(min(smooth.expr.ord.scaled), seq(-2, 2, length.out=99), max(smooth.expr.ord.scaled)),
    labels_row=rep("", nrow(smooth.expr.ord.scaled)),
    labels_col=rep("", ncol(smooth.expr.ord.scaled))
  )
}

#' Title
#'
#' @param time 
#' @param progression 
#'
#' @return
#' @export
#'
#' @examples
evaluate.trajectory <- function(time, progression) {
  stime <- sort(time)
  
  # if some timepoints contain multiple cells, randomly order these
  diff <- stime[-1]-stime[-length(stime)]
  min.diff <- min(diff[diff!=0])
  noises <- runif(length(time), 0, .01) * min.diff
  noised.time <- time + noises
  rank <- rank(noised.time)
  
  consistencies <- t(sapply(seq_along(rank), function(i) {
    left <- rank[[i]] > rank
    right <- rank[[i]] < rank
    up <- progression[[i]] > progression
    down <- progression[[i]] < progression
    
    leftup <- sum(left & up)
    leftdown <- sum(left & down)
    rightup <- sum(right & up)
    rightdown <- sum(right & down)
    
    one <- (leftup + rightdown)
    two <- (rightup + leftdown)
    
    consistency.direction1 <- one / (one + two) * 2 - 1
    consistency.direction2 <- two / (one + two) * 2 - 1
    
    c(consistency.direction1, consistency.direction2)
  }))
  
  consistency <- colMeans(consistencies)
  
  max(consistency)
}


#' Title
#'
#' @param space 
#' @param progression 
#' @param k 
#'
#' @return
#' @export
#' 
#' @importFrom class knn
#'
#' @examples
evaluate.space <- function(space, progression, k=5) {
  pred <- sapply(seq_len(nrow(space)), function(i) {
    as.integer(as.character(class::knn(space[-i,,drop=F], space[i,,drop=F], progression[-i], k=k)))
  })
  mean(progression==pred)
}

