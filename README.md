
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IsolationForest

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/JSzitas/IsolationForest.svg?branch=master)](https://travis-ci.org/JSzitas/IsolationForest)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/JSzitas/IsolationForest?branch=master&svg=true)](https://ci.appveyor.com/project/JSzitas/IsolationForest)
[![Codecov test
coverage](https://codecov.io/gh/JSzitas/IsolationForest/branch/master/graph/badge.svg)](https://codecov.io/gh/JSzitas/IsolationForest?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/IsolationForest)](https://CRAN.R-project.org/package=IsolationForest)
<!-- badges: end -->

IsolationForest implements Isolation Forests and Extended Isolation
Forests in R, optionally parallelized (for speed) using the **future**
framework.

Originally, this was cloned from Zelazny7/isofor
[GitHub](https://github.com/Zelazny7/isofor). This package, however,
also implements the Extended Isolation Forests of Hariri et al.(2019).
Further, there support for optional encoding of categorical variables
using the **categoryEncodings** package, and missing values can be
handled via Missingness Incorporated in Attributes splitting (Kwala et
al. (2008)).

**NOTE** The package does not explicitly handle factor splitting, as it
seems encoding factors might be a more reasonable approach for trees,
see ‘Sufficient Representations for Categorical Variables’ by
Johannemann et al. (2019) - the package currently supports these
encodings, hence the use of **categoryEncodings**.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JSzitas/IsolationForests")
```

<!-- Hopefully, once a CRAN release is made, the package will also be available via the usual route: -->

<!-- ``` r -->

<!-- install.packages("IsolationForests") -->

<!-- ``` -->

## Examples

Generate random data with anomalies:

``` r
set.seed(1071)
X <- rnorm( 500, 0, 1)
Y <- rnorm( 500, 0, 1)
replace_x <- sample(1:500, 20 )
replace_y <- sample(1:500, 50 )
X[replace_x] <- rnorm(20, mean = 3, sd = 2)
Y[replace_y] <- rnorm(50, mean = -4, sd = 1.5)

anomaly_indicator <- rep(0,500)
anomaly_indicator[replace_x] <- 1
anomaly_indicator[replace_y] <- 1
anomaly_indicator <- as.factor(anomaly_indicator)

test_data <- data.frame(X, Y, anomaly_indicator)


ggplot2::ggplot(data = test_data, ggplot2::aes( x = X,
                                                y = Y,
                                                colour = anomaly_indicator,
                                                shape = anomaly_indicator )) +
  ggplot2::geom_point(size = 1.9) +
  ggplot2::scale_colour_manual(name = "Anomaly", values = c("#2554C7","#E42217")) +
  ggplot2::scale_shape_manual(name = "Anomaly", values = c(15,17))
  
```

We have a total of

``` r
total_anomalous <- sum(unlist(anomaly_indicator == 1)) 
```

anomalous values. Keep that in mind for later.

Now try fitting an Isolation Forest:

``` r
library(IsolationForest)

fit <- isolationForest( X = test_data[,1:2], # we dont want column 3 here. 
                        n_trees = 1000,
                        Phi = 64, # subsampling rate for individual trees
                        parallel = TRUE, # defaults to future::plan("multiprocess")
                        future_plan = "multiprocess", # change this argument 
                                                      # to change the plan 
                        extension_level = 1, # how 'extended' should the trees be? 
                        vanilla = FALSE      # whether to fit an unextended, original 
                                             # isolation forest
                        )
```

Then to get the anomaly scores we just call

``` r
 scored_data <- predict.isolationForest(fit, test_data[,1:2])
 
```

We can additionaly generate 2 dimensional contour plots by calling

``` r
anomaly_plot( x = "X",
              y = "Y",
              forest = fit,
              data = test_data[,1:2], 
              contour = TRUE )
```

Or we can plot the individual point, classified as anomalous or not

``` r
anomaly_plot( x = "X",
              y = "Y",
              forest = fit,
              data = test_data[,1:2], 
              contour = FALSE,
              contamination = 0.15 # we have contaminated total_anomalous/nrow(test_data),
                                   # observations 
              )
```
