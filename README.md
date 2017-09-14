<!-- github markdown built using 
rmarkdown::render("README.Rmd", output_format = "md_document")
-->
SCORPIUS
========

[![Build Status](https://travis-ci.org/rcannood/SCORPIUS.png?branch=master)](https://travis-ci.org/rcannood/SCORPIUS)

**SCORPIUS an unsupervised approach for inferring developmental chronologies from single-cell RNA sequencing data.** In comparison to similar approaches, it has three main advantages:

-   **It accurately reconstructs linear dynamic processes.** The performance was evaluated using a quantitative evaluation pipeline and ten single-cell RNA sequencing datasets.

-   **It automatically identifies marker genes, speeding up knowledge discovery.**

-   **It is fully unsupervised.** Prior knowledge of the relevant marker genes or cellular states of individual cells is not required.

News:

-   A preprint is available on [bioRxiv](http://biorxiv.org/content/early/2016/10/07/079509).

-   Check out our [review](http://onlinelibrary.wiley.com/doi/10.1002/eji.201646347/full) on Trajectory Inference methods!

Installing SCORPIUS
-------------------

You can install the latest version from github with

``` r
devtools::install_github("rcannood/SCORPIUS", build_vignettes=TRUE)
```

<!--
You can install:

* the latest released version from CRAN with

    ```R
    install.packages("SCORPIUS")
    ```

* the latest development version from github with

    ```R
    devtools::install_github("rcannood/SCORPIUS", build_vignettes=TRUE)
    ```
-->
If you encounter a bug, please file a minimal reproducible example on the [issues](https://github.com/rcannood/SCORPIUS/issues) page.

Learning SCORPIUS
-----------------

To get started, read the introductory example below, or read one of the vignettes containing more elaborate examples:

-   [Investigating differentiating dendritic cell progenitors](vignettes/ginhoux.md): `vignette("ginhoux", package="SCORPIUS")`
-   [Inferring trajectories from simulated data](vignettes/simulated-data.md): `vignette("simulated-data", package="SCORPIUS")`

Introductory example
--------------------

This section describes the main workflow of SCORPIUS without going in depth in the R code. For a more detailed explanation, see the vignettes listed below.

To start using SCORPIUS, simply write:

``` r
library(SCORPIUS)
```

The `ginhoux` dataset (See Schlitzer et al. 2015) contains 248 dendritic cell progenitors in one of three cellular cellular states: MDP, CDP or PreDC.

``` r
data(ginhoux)
expression <- ginhoux$expression
group_name <- ginhoux$sample_info$group_name
```

With the following code, SCORPIUS reduces the dimensionality of the dataset and provides a visual overview of the dataset. In this plot, cells that are similar in terms of expression values will be placed closer together than cells with dissimilar expression values.

``` r
dist <- correlation_distance(expression)

# filter outliers
filt <- outlier_filter(dist)
expression <- expression[filt, ]
group_name <- group_name[filt]
dist <- dist[filt, filt]

# reduce dimensionality
space <- reduce_dimensionality(dist)
draw_trajectory_plot(space, group_name)
```

![](README_files/figure-markdown_github/reduce%20dimensionality-1.png)

To infer and visualise a trajectory through these cells, run:

``` r
traj <- infer_trajectory(space)
draw_trajectory_plot(space, group_name, traj$path)
```

![](README_files/figure-markdown_github/infer%20trajectory-1.png)

To identify candidate marker genes

``` r
# warning: setting num_permutations to 10 requires a long time (~30min) to run!
# set it to 0 and define a manual cutoff for the genes (e.g. top 200) for a much shorter execution time.
gimp <- gene_importances(
  expression, 
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 
```

Select the most important genes, scale its expession, cluster them into modules, and visualise it.

``` r
gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(expression[,gene_sel])
modules <- extract_modules(expr_sel, traj$time)
```


``` r
# data is already quantile scaled
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, scale_features = F)
```

![](README_files/figure-markdown_github/visualise%20tafs-1.png)

Related approaches
------------------

-   [Wanderlust](http://www.c2b2.columbia.edu/danapeerlab/html/wanderlust.html)
-   [Monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html)
-   [Waterfall](http://dx.doi.org/10.1016/j.stem.2015.07.013)
-   [Embeddr](https://github.com/kieranrcampbell/embeddr)

References
----------

Schlitzer, Andreas, V Sivakamasundari, Jinmiao Chen, Hermi Rizal Bin Sumatoh, Jaring Schreuder, Josephine Lum, Benoit Malleret, et al. 2015. “Identification of cDC1- and cDC2-committed DC progenitors reveals early lineage priming at the common DC progenitor stage in the bone marrow.” *Nature Immunology* 16 (7): 718–26. doi:[10.1038/ni.3200](https://doi.org/10.1038/ni.3200).
