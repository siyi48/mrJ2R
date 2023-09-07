
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Codes Description

<!-- badges: start -->
<!-- badges: end -->

The goal of mrJ2R is to obtain the proposed estimators to evaluate the
average treatment effect (ATE) under jump-to-reference (J2R) based on
the paper **Multiply robust estimators in longitudinal studies with
missing data under control-based imputation**.

## Installation

You can install the development version of rpsftmPDT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("siyi48/mrJ2R")
#> Skipping install of 'mrJ2R' from a github remote, the SHA1 (1981fb41) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Usage

In the cross-sectional study, the main function `jtr.1time` provides the
ATE estimator

``` r
library(mrJ2R)
data <- dat1
formula.om <- y ~ z1 + z2 + z3 + z4 + z5 + a +
  a:(z1 + z2 + z3 + z4 + z5)
formula.ps <- a ~ z1 + z2 + z3 + z4 + z5
formula.rp <- r ~ z1 + z2 + z3 + z4 + z5 + a +
  a:(z1 + z2 + z3 + z4 + z5)
mat.cal <- data.matrix(data[,paste0("z", 1:5)]) # calibration matrix
res.1time <- jtr.1time(formula.ps = formula.ps,
                       formula.rp = formula.rp,
                       formula.om = formula.om,
                       data = data, mat.cal = mat.cal,
                       type = c("tr", "tr.cal", "rpom"))
res.1time
#> $est
#>           tr       tr.cal         rpom 
#> 0.0001152698 0.0070087678 0.0089585440 
#> 
#> $ve
#>          tr      tr.cal 
#> 0.003680655 0.003028245
```
