``` r
library(SCORPIUS)
library(dplyr)
library(tidyr)
```

Loading in the data
-------------------

Load in the dataset from Schlitzer et al.

``` r
data(ginhoux.dataset)
```

This matrix contains the expression data

``` r
counts <- ginhoux.dataset$counts
dim(counts)
```

    ## [1]   248 15752

This is info pertaining the cells

``` r
sample.info <- ginhoux.dataset$sample.info
head(sample.info)
```

    ##            population
    ## SRR1558744        CDP
    ## SRR1558745        CDP
    ## SRR1558746        CDP
    ## SRR1558747        CDP
    ## SRR1558748        CDP
    ## SRR1558749        CDP

``` r
celltypes <- ginhoux.dataset$celltypes %>% mutate(i=seq_along(parent)) %>% select(-parent)
celltypes
```

    ##   population i
    ## 1        MDP 1
    ## 2        CDP 2
    ## 3      PreDC 3

We construct a vector of the progression of cells as integers and as factors

``` r
progression <- sample.info %>% left_join(celltypes, by=colnames(celltypes)[-ncol(celltypes)]) %>% .$i

levels <- celltypes %>% select(-i) %>% apply(1, paste, collapse="_")
progression.str <- factor(sample.info[,colnames(celltypes)[-ncol(celltypes)],drop=F] %>% apply(1, paste, collapse="_"), levels)

names(progression.str) <- names(progression) <- rownames(counts)
```

SCORPIUS
--------

Calculate distance

``` r
dist <- calculate.distance(counts)
```

Calculate outliers

``` r
out <- calculate.outlier(dist)
filter <- !out$is.outlier
```

Filter away outliers

``` r
counts <- counts[filter,,drop=F]
dist <- dist[filter, filter]
sample.info <- sample.info[filter,,drop=F]
progression.str <- progression.str[filter]
progression <- progression[filter]
```

Reduce dimensionality

``` r
space <- reduce.dimensionality(dist, ndim=3)
evaluate.space(space, progression)
```

    ## [1] 0.9628099

Different plotting options

``` r
plot.dimensionality.reduction(space, colour = progression.str)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
plot.dimensionality.reduction(space, colour = progression.str, contour=T)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-2.png)

``` r
plot.dimensionality.reduction(space, colour = progression.str, contour=T)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-3.png)

``` r
plot.dimensionality.reduction(space, contour=T)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-4.png)

Infer trajectory

``` r
trajectory <- infer.trajectory(space, k=4)
evaluate.trajectory(trajectory$time, progression)
```

    ## [1] 0.9589824

Different plotting options

``` r
plot.trajectory(space, trajectory$initial.path, colour=progression.str, contour=T)
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
plot.trajectory(space, trajectory$final.path, colour=progression.str, contour=T)
```

![](README_files/figure-markdown_github/unnamed-chunk-12-2.png)

``` r
plot.trajectory.density(trajectory$time, progression.str)
```

![](README_files/figure-markdown_github/unnamed-chunk-12-3.png)

Reverse the trajectory if it is the wrong way around

``` r
prog.means <- tapply(trajectory$time, progression.str, mean)
if (prog.means[[last(levels)]] < prog.means[[first(levels)]]) {
  trajectory <- reverse.trajectory(trajectory)
}
time <- trajectory$time
plot.trajectory.density(trajectory$time, progression.str)
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

Find trajectory aligned genes (TAGs)

``` r
tags <- find.trajectory.aligned.genes(counts, trajectory$time, q.value.cutoff = 1e-10, mc.cores = 8)
head(tags$p.values, 20)
```

    ##             gene      p.value      q.value is.tag           category
    ## 1            Mpo 5.477392e-82 8.627988e-78   TRUE         p <= 1e-40
    ## 2          Prtn3 7.036260e-72 5.541758e-68   TRUE         p <= 1e-40
    ## 3          H2-Aa 1.483380e-59 7.788732e-56   TRUE         p <= 1e-40
    ## 4           Cd74 1.288486e-55 5.074059e-52   TRUE         p <= 1e-40
    ## 5            Cd7 3.088776e-49 6.950629e-46   TRUE         p <= 1e-40
    ## 6           Ctsg 2.790387e-49 6.950629e-46   TRUE         p <= 1e-40
    ## 7         H2-Ab1 2.927322e-49 6.950629e-46   TRUE         p <= 1e-40
    ## 8           Cd34 2.174884e-44 4.282347e-41   TRUE         p <= 1e-40
    ## 9         H2-Eb1 2.672143e-43 4.676844e-40   TRUE 1e-40 < p <= 1e-20
    ## 10          Ly6d 2.843103e-42 4.478455e-39   TRUE 1e-40 < p <= 1e-20
    ## 11      AK152437 4.747042e-38 6.231284e-35   TRUE 1e-40 < p <= 1e-20
    ## 12          Cd93 4.516503e-38 6.231284e-35   TRUE 1e-40 < p <= 1e-20
    ## 13 2810417H13Rik 1.412088e-37 1.711017e-34   TRUE 1e-40 < p <= 1e-20
    ## 14         Top2a 3.300674e-37 3.713730e-34   TRUE 1e-40 < p <= 1e-20
    ## 15        Igfbp4 8.279744e-36 8.694835e-33   TRUE 1e-40 < p <= 1e-20
    ## 16          Cst7 7.813861e-35 7.692747e-32   TRUE 1e-40 < p <= 1e-20
    ## 17       Siglech 1.492379e-34 1.382821e-31   TRUE 1e-40 < p <= 1e-20
    ## 18         Plbd1 5.088580e-34 4.453073e-31   TRUE 1e-40 < p <= 1e-20
    ## 19         Stmn1 2.631893e-33 2.181978e-30   TRUE 1e-40 < p <= 1e-20
    ## 20         Tubb5 7.878000e-33 6.204713e-30   TRUE 1e-40 < p <= 1e-20

Group genes into modules

``` r
modules <- find.modules(tags$smooth.expression, tags$genes)
```

    ##  ..done.
    ##  mergeCloseModules: Merging modules whose distance is less than 0.3
    ##    Calculating new MEs...

``` r
plot.modules.heatmap(tags$smooth.expression, tags$genes, progression.str, time, modules)
```

![](README_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
plot.modules.heatmap(counts, tags$genes, progression.str, time, modules)
```

![](README_files/figure-markdown_github/unnamed-chunk-15-2.png)
