<!-- built using 
render("simulated-data.Rmd", output_format = "all") 
-->
In this vignette, SCORPIUS is used to infer a trajectory through cells in artificial single-cell RNA-seq data. Note that the dataset is generated in a very naive manner and is only meant to be used for demonstration purposes, not for evaluating trajectory inference methods.

Simulate expression data
------------------------

Expression values for 384 cells and 500 genes is generated as follows.

``` r
library(SCORPIUS)
dataset <- generate.dataset(type="poly", num.genes=500, num.samples=384, num.groups=4)
```

The resulting dataset is a list containing a matrix named `expression` and a data frame named `sample.info`.

`expression` is a 384-by-500 matrix containing the expression values of all the cells and all the genes.

``` r
dataset$expression[1:6, 1:6]
```

    ##            Gene1    Gene2    Gene3     Gene4     Gene5    Gene6
    ## Sample1 1.177482 0.000000 8.369938  3.198941  9.533305 4.754622
    ## Sample2 6.099312 0.000000 9.635797  5.251118  5.632712 5.571352
    ## Sample3 0.000000 5.771299 4.991685  8.092176  5.281727 8.220894
    ## Sample4 0.000000 3.530891 7.490812  4.176563  2.213085 6.726484
    ## Sample5 0.000000 0.000000 6.880374  6.221958  0.000000 5.976617
    ## Sample6 2.545596 3.495439 9.126022 11.032961 12.759797 5.760015

`sample.info` is a data frame with the metadata of the cells, containing only the group each cell belongs to.

``` r
head(dataset$sample.info)
```

    ##         group.name
    ## Sample1    Group 1
    ## Sample2    Group 1
    ## Sample3    Group 1
    ## Sample4    Group 1
    ## Sample5    Group 1
    ## Sample6    Group 1

In order to infer a trajectory through this data, SCORPIUS first reduces the dimensionality of the dataset.

Reduce dimensionality of the dataset
------------------------------------

SCORPIUS uses classical Torgerson multi-dimensional scaling to reduce the dataset to three dimensions. In short, this technique attempts to place the cells in a space such that the distance between any two points in that space approximates the original distance between the two cells as well as possible.

The distance between any two samples is defined as their correlation distance, namely `1 - (cor(x, y)+1)/2`. The distance matrix is calculated as follows:

``` r
dist <- correlation.distance(dataset$expression)
```

`dist` is a 384-by-384 matrix, with values ranging from 0 to 1.

``` r
dim(dist)
```

    ## [1] 384 384

``` r
plot(density(dist))
```

![](simulated-data_files/figure-markdown_github/unnamed-chunk-5-1.png)

The reduced space is constructed as follows:

``` r
space <- reduce.dimensionality(dist, ndim=3)
```

The new space is a 384-by-3 matrix, and can be visualised as follows:

``` r
draw.trajectory.plot(space)
```

![](simulated-data_files/figure-markdown_github/unnamed-chunk-7-1.png)

Looking at this plot, it seems that the cells in this dataset are involved in a dynamic process.

In addition, if a property of the cells (e.g. cell type) is known, it can be used to colour the plot using the `progression.group` parameter.

In this case, the underlying groups of each cell were also given:

``` r
group <- dataset$sample.info$group.name
draw.trajectory.plot(space, progression.group = group)
```

![](simulated-data_files/figure-markdown_github/unnamed-chunk-8-1.png)

Inferring a trajectory through the cells
----------------------------------------

The main goal of SCORPIUS is to infer a trajectory through the cells, and orden the cells according to the inferred timeline.

SCORPIUS infers a trajectory through several intermediate steps, which are all executed as follows:

``` r
traj <- infer.trajectory(space, k=4)
```

The result is a list containing the final trajectory `final.path` and the inferred timeline for each sample `time`.

The trajectoryy can be visualised with respect to the samples by passing it to `draw.trajectory.plot`:

``` r
draw.trajectory.plot(space, progression.group = group, path = traj$final.path)
```

![](simulated-data_files/figure-markdown_github/unnamed-chunk-10-1.png)

Finding genes of interest
-------------------------

Coming soon...
