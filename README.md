<!-- README.md is generated from README.Rmd. Please edit that file -->
SCORPIUS
========

[![Build Status](https://travis-ci.org/rcannood/SCORPIUS.png?branch=master)](https://travis-ci.org/rcannood/SCORPIUS)

SCORPIUS an unsupervised approach for inferring developmental chronologies from single-cell RNA sequencing data. In comparison to similar approaches, it has three main advantages:

-   It accurately reconstructs trajectories for a wide variety of dynamic cellular processes. The performance was evaluated using a new, quantitative evaluation pipeline, comparing the performance of current state-of-the-art techniques on 10 publicly available single-cell RNA sequencing datasets.

-   It identifies marker genes. By automatically identifying possible marker genes relevant for the dynamic process under investigation, SCORPIUS speeds up

-   It is fully unsupervised. Prior knowledge of the relevant marker genes or cellular states of individual cells is not required, and thus the conclusions drawn from the SCORPIUS output are unbiased.

You can install the latest version from github with

``` r
devtools::install_github("rcannood/SCORPIUS")
```

<!--
You can install:

* the latest released version from CRAN with

    ```R
    install.packages("SCORPIUS")
    ```

* the latest development version from github with

    ```R
    devtools::install_github("rcannood/SCORPIUS")
    ```
-->
If you encounter a clear bug, please file a minimal reproducible example on [github](https://github.com/rcannood/SCORPIUS/issues).

Learning SCORPIUS
-----------------

To get started, read the notes below or read the ~~[intro vignette](): `vignette("introduction", package="SCORPIUS")`~~ (Under construction).

More vignettes with more elaborate examples are also available:

-   ~~[Investigating differentiating dendritic cell progenitors](): `vignette("ginhoux", package="SCORPIUS")`~~ (Under construction)
-   [Inferring trajectories from simulated data](vignettes/simulated-data.md): `vignette("simulated-data", package="SCORPIUS")`

Minimal example
---------------

To start using SCORPIUS, simply write:

``` r
library(SCORPIUS)
```

The `ginhoux` dataset (See Schlitzer et al. 2015) contains 248 cell progenitors in one of three cellular cellular states: MDP, CDP or PreDC.

``` r
data(ginhoux)
```

With the following code, SCORPIUS reduces the dimensionality of the dataset and provides a visual overview of the dataset. In this plot, cells that are similar in terms of expression values will be placed closer together than cells with dissimilar expression values.

``` r
dist <- correlation.distance(ginhoux$expression)
space <- reduce.dimensionality(dist)
group.name <- ginhoux$sample.info$group.name
draw.trajectory.plot(space, group.name, contour = T)
```

![](README_files/figure-markdown_github/reduce%20dimensionality-1.png)
 To infer and visualise a trajectory through these cells, run:

``` r
traj <- infer.trajectory(space)
draw.trajectory.plot(space, group.name, traj$final.path, contour = T)
```

![](README_files/figure-markdown_github/infer%20trajectory-1.png)
 Finally, to identify and visualise candidate marker genes, execute the following code:

``` r
tafs <- find.trajectory.aligned.features(ginhoux$expression, traj$time)
expr.tafs <- tafs$smooth.x[,tafs$tafs]
modules <- extract.modules(expr.tafs)
draw.trajectory.heatmap(expr.tafs, traj$time, group.name, modules)
```

![](README_files/figure-markdown_github/find%20tafs-1.png)

Related approaches
------------------

-   [Monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html)
-   [Waterfall](http://dx.doi.org/10.1016/j.stem.2015.07.013)
-   [Embeddr](https://github.com/kieranrcampbell/embeddr)

Schlitzer, Andreas, V Sivakamasundari, Jinmiao Chen, Hermi Rizal Bin Sumatoh, Jaring Schreuder, Josephine Lum, Benoit Malleret, et al. 2015. “Identification of cDC1- and cDC2-committed DC progenitors reveals early lineage priming at the common DC progenitor stage in the bone marrow.” *Nature Immunology* 16 (7): 718–26. doi:[10.1038/ni.3200](http://dx.doi.org/10.1038/ni.3200).
