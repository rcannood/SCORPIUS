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

-   ~~[Investigating differentiating dendritic cells](): `vignette("ginhoux", package="SCORPIUS")`~~ (Under construction)
-   [Inferring trajectories from simulated data](https://github.com/rcannood/SCORPIUS/blob/master/vignettes/simulated-data.md): `vignette("simulated-data", package="SCORPIUS")`

<!--
## Learning SCORPIUS

To get started, read the notes below, then read the intro vignette: `vignette("introduction", package = "dplyr")`. To make the most of dplyr, I also recommend that you familiarise yourself with the principles of [tidy data](http://vita.had.co.nz/papers/tidy-data.html): this will help you get your data into a form that works well with dplyr, ggplot2 and R's many modelling functions.

If you need more, help I recommend the following (paid) resources:

* [dplyr](https://www.datacamp.com/courses/dplyr) on datacamp, by Garrett
  Grolemund. Learn the basics of dplyr at your own pace in this interactive 
  online course.
  
* [Introduction to Data Science with R](http://shop.oreilly.com/product/0636920034834.do): 
  How to Manipulate, Visualize, and Model Data with the R Language, by Garrett
  Grolemund. This O'Reilly video series will teach you the basics needed to be
  an effective analyst in R.

## Key data structures

The key object in dplyr is a _tbl_, a representation of a tabular data structure.
Currently `dplyr` supports:

* data frames
* [data tables](https://github.com/Rdatatable/data.table/wiki)
* [SQLite](http://sqlite.org/)
* [PostgreSQL](http://www.postgresql.org/)/[Redshift](http://aws.amazon.com/redshift/)
* [MySQL](http://www.mysql.com/)/[MariaDB](https://mariadb.com/)
* [Bigquery](https://developers.google.com/bigquery/)
* [MonetDB](http://www.monetdb.org/)
* data cubes with arrays (partial implementation)

You can create them as follows:

```R
{r, message = FALSE}
library(dplyr) # for functions
library(nycflights13) # for data
flights

# Caches data in local SQLite db
flights_db1 <- tbl(nycflights13_sqlite(), "flights")

# Caches data in local postgres db
flights_db2 <- tbl(nycflights13_postgres(), "flights")
```

Each tbl also comes in a grouped variant which allows you to easily perform operations "by group":

```R
{r}
carriers_df  <- flights %>% group_by(carrier)
carriers_db1 <- flights_db1 %>% group_by(carrier)
carriers_db2 <- flights_db2 %>% group_by(carrier)
```

## Single table verbs

`dplyr` implements the following verbs useful for data manipulation:

* `select()`: focus on a subset of variables
* `filter()`: focus on a subset of rows
* `mutate()`: add new columns
* `summarise()`: reduce each group to a smaller number of summary statistics
* `arrange()`: re-order the rows

They all work as similarly as possible across the range of data sources. The main difference is performance:

```R
{r}
system.time(carriers_df %>% summarise(delay = mean(arr_delay)))
system.time(carriers_db1 %>% summarise(delay = mean(arr_delay)) %>% collect())
system.time(carriers_db2 %>% summarise(delay = mean(arr_delay)) %>% collect())
```

Data frame methods are much much faster than the plyr equivalent. The database methods are slower, but can work with data that don't fit in memory.

```R
{r}
system.time(plyr::ddply(flights, "carrier", plyr::summarise,
  delay = mean(arr_delay, na.rm = TRUE)))
```

### `do()`

As well as the specialised operations described above, `dplyr` also provides the generic `do()` function which applies any R function to each group of the data.

Let's take the batting database from the built-in Lahman database. We'll group it by year, and then fit a model to explore the relationship between their number of at bats and runs:

```R
{r}
by_year <- lahman_df() %>% 
  tbl("Batting") %>%
  group_by(yearID)
by_year %>% 
  do(mod = lm(R ~ AB, data = .))
```

Note that if you are fitting lots of linear models, it's a good idea to use `biglm` because it creates model objects that are considerably smaller:

```R
{r}
by_year %>% 
  do(mod = lm(R ~ AB, data = .)) %>%
  object.size() %>%
  print(unit = "MB")

by_year %>% 
  do(mod = biglm::biglm(R ~ AB, data = .)) %>%
  object.size() %>%
  print(unit = "MB")
```

## Multiple table verbs

As well as verbs that work on a single tbl, there are also a set of useful verbs that work with two tbls at a time: joins and set operations.

dplyr implements the four most useful joins from SQL:

* `inner_join(x, y)`: matching x + y
* `left_join(x, y)`: all x + matching y
* `semi_join(x, y)`: all x with match in y
* `anti_join(x, y)`: all x without match in y

And provides methods for:

* `intersect(x, y)`: all rows in both x and y
* `union(x, y)`: rows in either x or y
* `setdiff(x, y)`: rows in x, but not y

## Plyr compatibility

You'll need to be a little careful if you load both plyr and dplyr at the same time. I'd recommend loading plyr first, then dplyr, so that the faster dplyr functions come first in the search path. By and large, any function provided by both dplyr and plyr works in a similar way, although dplyr functions tend to be faster and more general.

## Related approaches

* [Blaze](http://blaze.pydata.org)
* [|Stat](http://oldwww.acm.org/perlman/stat/)
* [Pig](http://dx.doi.org/10.1145/1376616.1376726)
-->
