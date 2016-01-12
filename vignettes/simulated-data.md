<!-- built using 
render("vignettes/simulated-data.Rmd", output_format = "all") 
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

    ##             Gene1    Gene2    Gene3    Gene4     Gene5     Gene6
    ## Sample1 5.9320775 0.000000 0.000000 1.931591  8.041771  3.221728
    ## Sample2 0.0000000 0.000000 0.000000 7.999034  5.781027  1.895666
    ## Sample3 0.2835773 0.000000 4.882344 2.410108  7.248168  9.296099
    ## Sample4 3.9561486 5.132754 0.000000 8.795270 10.075464  6.100947
    ## Sample5 3.7062338 0.000000 0.000000 0.000000  8.345310  5.400626
    ## Sample6 0.0000000 7.079570 0.000000 5.102890  0.000000 11.879792

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
group.name <- dataset$sample.info$group.name
draw.trajectory.plot(space, progression.group = group.name)
```

![](simulated-data_files/figure-markdown_github/unnamed-chunk-8-1.png)

Inferring a trajectory through the cells
----------------------------------------

The main goal of SCORPIUS is to infer a trajectory through the cells, and orden the cells according to the inferred timeline.

SCORPIUS infers a trajectory through several intermediate steps, which are all executed as follows:

``` r
traj <- infer.trajectory(space)
```

The result is a list containing the final trajectory `final.path` and the inferred timeline for each sample `time`.

The trajectory can be visualised with respect to the samples by passing it to `draw.trajectory.plot`:

``` r
draw.trajectory.plot(space, progression.group = group.name, path = traj$final.path)
```

![](simulated-data_files/figure-markdown_github/unnamed-chunk-10-1.png)

Finding candidate marker genes
------------------------------

Finally, to identify and visualise candidate marker genes, execute the following code:

``` r
tafs <- find.trajectory.aligned.features(dataset$expression, traj$time)
expr.tafs <- tafs$smooth.x[,tafs$tafs]
modules <- extract.modules(expr.tafs)
draw.trajectory.heatmap(expr.tafs, traj$time, group.name, modules)
```

![](simulated-data_files/figure-markdown_github/find%20tafs-1.png)
