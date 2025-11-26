# scRNA-seq data of dendritic cell progenitors.

This dataset contains the expression values of the top 2000 most
variable genes for 248 dendritic cell progenitors. Each cell is in one
of three maturation stages: MDP, CDP or PreDC. The levels of the factor
in `sample.info` are ordered according to the maturation process.

The number of genes had to be reduced specifically for reducing the
package size of SCORPIUS. Use the following code to download the
original data:

    download.file("https://github.com/rcannood/SCORPIUS/raw/master/data-raw/ginhoux_orig.rds", destfile = "local.rds")
    ginhoux <- readRDS("local.rds")
    # do something with ginhoux

## Usage

``` r
ginhoux
```

## Format

A list containing two data frames, `expression` (248x2000) and
`sample_info` (248x1).

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60783>

## References

Schlitzer A, Sivakamasundari V, Chen J, Sumatoh HR et al. Identification
of cDC1- and cDC2-committed DC progenitors reveals early lineage priming
at the common DC progenitor stage in the bone marrow. Nat Immunol 2015
Jul;16(7):718-28. PMID: 26054720

## See also

[SCORPIUS](rcannood.github.io/SCORPIUS/reference/SCORPIUS-package.md)
