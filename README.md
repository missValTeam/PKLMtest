# PKLMtest: testing MCAR with classification


## Overview

PKLMtest is a package intended to provide a framework for testing MCAR with classification. It implements the test  described in [Michel, Naef, Spohn and Meinshausen. 2021](https://arxiv.org/abs/2109.10150) . Examples of use of the library are shown below.


## Installation

To install the package from github you can run

``` r
install.packages("devtools")
devtools::install_github("missValTeam/PKLMtest")
```

## Examples: 

```r
n <- 500 
X <- cbind(rnorm(n),rnorm(n))
X.NA <- X
X.NA[,1] <- ifelse(stats::runif(n)<=0.2, NA, X[,1])
pval <- PKLMtest(X)
```


## Issues

To report an issue, please use the [issue tracker](https://github.com/missValTeam/PKLMtest/issues) on github.com.
