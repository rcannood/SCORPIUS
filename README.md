``` r
library(SCORPIUS)
```

    ## Loading required package: DBI

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
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
progression
```

    ##   [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ##  [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ##  [71] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3
    ## [106] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    ## [141] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    ## [176] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [211] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [246] 1 1 1

``` r
levels <- celltypes %>% select(-i) %>% apply(1, paste, collapse="_")
progression.str <- factor(sample.info[,colnames(celltypes)[-ncol(celltypes)],drop=F] %>% apply(1, paste, collapse="_"), levels)
progression.str
```

    ## SRR1558744 SRR1558745 SRR1558746 SRR1558747 SRR1558748 SRR1558749 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558750 SRR1558751 SRR1558752 SRR1558753 SRR1558754 SRR1558755 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558756 SRR1558757 SRR1558758 SRR1558759 SRR1558760 SRR1558761 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558762 SRR1558763 SRR1558764 SRR1558765 SRR1558766 SRR1558767 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558768 SRR1558769 SRR1558770 SRR1558771 SRR1558772 SRR1558773 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558774 SRR1558775 SRR1558776 SRR1558777 SRR1558778 SRR1558779 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558780 SRR1558781 SRR1558782 SRR1558783 SRR1558784 SRR1558785 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558786 SRR1558787 SRR1558788 SRR1558790 SRR1558791 SRR1558792 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558793 SRR1558794 SRR1558795 SRR1558796 SRR1558797 SRR1558798 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558799 SRR1558800 SRR1558801 SRR1558802 SRR1558803 SRR1558804 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558805 SRR1558806 SRR1558807 SRR1558808 SRR1558809 SRR1558810 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558811 SRR1558812 SRR1558813 SRR1558814 SRR1558815 SRR1558816 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558817 SRR1558818 SRR1558819 SRR1558820 SRR1558821 SRR1558822 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558823 SRR1558824 SRR1558825 SRR1558826 SRR1558827 SRR1558828 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558829 SRR1558830 SRR1558831 SRR1558832 SRR1558833 SRR1558834 
    ##        CDP        CDP        CDP        CDP        CDP        CDP 
    ## SRR1558835 SRR1558836 SRR1558837 SRR1558838 SRR1558839 SRR1558840 
    ##        CDP        CDP        CDP        CDP        CDP      PreDC 
    ## SRR1558841 SRR1558842 SRR1558843 SRR1558844 SRR1558845 SRR1558846 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558847 SRR1558848 SRR1558849 SRR1558850 SRR1558851 SRR1558852 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558853 SRR1558854 SRR1558855 SRR1558856 SRR1558857 SRR1558858 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558859 SRR1558860 SRR1558861 SRR1558862 SRR1558863 SRR1558864 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558865 SRR1558866 SRR1558867 SRR1558868 SRR1558869 SRR1558870 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558871 SRR1558872 SRR1558873 SRR1558874 SRR1558875 SRR1558876 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558877 SRR1558878 SRR1558879 SRR1558880 SRR1558881 SRR1558882 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558883 SRR1558884 SRR1558885 SRR1558886 SRR1558887 SRR1558888 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558889 SRR1558890 SRR1558891 SRR1558892 SRR1558893 SRR1558894 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558895 SRR1558896 SRR1558897 SRR1558898 SRR1558899 SRR1558900 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558901 SRR1558902 SRR1558903 SRR1558904 SRR1558905 SRR1558906 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558907 SRR1558908 SRR1558909 SRR1558910 SRR1558911 SRR1558912 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558913 SRR1558914 SRR1558915 SRR1558916 SRR1558917 SRR1558918 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558919 SRR1558920 SRR1558921 SRR1558922 SRR1558923 SRR1558924 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558925 SRR1558926 SRR1558927 SRR1558928 SRR1558929 SRR1558930 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC      PreDC 
    ## SRR1558931 SRR1558932 SRR1558933 SRR1558934 SRR1558935 SRR1558936 
    ##      PreDC      PreDC      PreDC      PreDC      PreDC        MDP 
    ## SRR1558937 SRR1558938 SRR1558939 SRR1558940 SRR1558941 SRR1558942 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558943 SRR1558944 SRR1558945 SRR1558946 SRR1558947 SRR1558948 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558949 SRR1558950 SRR1558951 SRR1558952 SRR1558953 SRR1558954 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558955 SRR1558956 SRR1558957 SRR1558958 SRR1558959 SRR1558960 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558961 SRR1558962 SRR1558963 SRR1558964 SRR1558965 SRR1558966 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558967 SRR1558968 SRR1558969 SRR1558970 SRR1558971 SRR1558972 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558973 SRR1558974 SRR1558975 SRR1558976 SRR1558977 SRR1558978 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558979 SRR1558980 SRR1558981 SRR1558983 SRR1558984 SRR1558985 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558987 SRR1558988 SRR1558989 SRR1558990 SRR1558991 SRR1558992 
    ##        MDP        MDP        MDP        MDP        MDP        MDP 
    ## SRR1558993 SRR1558994 
    ##        MDP        MDP 
    ## Levels: MDP CDP PreDC

``` r
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
trajectory <- infer.trajectory(space, k=5)
evaluate.trajectory(trajectory$time, progression)
```

    ## [1] 0.9573506

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
    ## 1            Mpo 1.309449e-84 2.062644e-80   TRUE         p <= 1e-40
    ## 2           Cd74 1.044825e-58 5.486026e-55   TRUE         p <= 1e-40
    ## 3          H2-Aa 7.651678e-59 5.486026e-55   TRUE         p <= 1e-40
    ## 4          Prtn3 7.242680e-57 2.852167e-53   TRUE         p <= 1e-40
    ## 5         H2-Ab1 8.936230e-51 2.815270e-47   TRUE         p <= 1e-40
    ## 6           Ctsg 3.502673e-49 9.195683e-46   TRUE         p <= 1e-40
    ## 7         H2-Eb1 3.986941e-47 8.971757e-44   TRUE         p <= 1e-40
    ## 8            Cd7 2.642793e-44 5.203659e-41   TRUE         p <= 1e-40
    ## 9       AK152437 1.203985e-41 2.107242e-38   TRUE 1e-40 < p <= 1e-20
    ## 10          Cd93 4.207735e-41 6.628024e-38   TRUE 1e-40 < p <= 1e-20
    ## 11          Cd34 2.515557e-40 3.602278e-37   TRUE 1e-40 < p <= 1e-20
    ## 12         Plbd1 3.111277e-38 4.084070e-35   TRUE 1e-40 < p <= 1e-20
    ## 13        Igfbp4 3.050987e-36 3.696858e-33   TRUE 1e-40 < p <= 1e-20
    ## 14          Cst7 1.760455e-35 1.980764e-32   TRUE 1e-40 < p <= 1e-20
    ## 15         Tubb5 3.168550e-35 3.327400e-32   TRUE 1e-40 < p <= 1e-20
    ## 16         Top2a 1.110690e-34 1.093474e-31   TRUE 1e-40 < p <= 1e-20
    ## 17       Siglech 2.215512e-34 2.052867e-31   TRUE 1e-40 < p <= 1e-20
    ## 18          Actb 5.011950e-34 4.386014e-31   TRUE 1e-40 < p <= 1e-20
    ## 19 2810417H13Rik 1.878208e-33 1.557133e-30   TRUE 1e-40 < p <= 1e-20
    ## 20           Myb 2.824244e-33 2.224375e-30   TRUE 1e-40 < p <= 1e-20

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
