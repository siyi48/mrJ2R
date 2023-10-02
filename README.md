
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
#> Downloading GitHub repo siyi48/mrJ2R@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/9h/trht1flx5gj0n8ntb4n7jlr40000gn/T/RtmpZd78I8/remotes1f2759e1441d/siyi48-mrJ2R-5f43eb7/DESCRIPTION’ ... OK
#> * preparing ‘mrJ2R’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘mrJ2R_0.0.0.9000.tar.gz’
#> Installing package into '/private/var/folders/9h/trht1flx5gj0n8ntb4n7jlr40000gn/T/Rtmp1zjFLn/temp_libpath199a63b161de'
#> (as 'lib' is unspecified)
```

## Usage

In the cross-sectional study, the main function `jtr.1time` provides the
ATE estimator.

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

In the longitudinal study, the main function `jtr.longi` provides the
ATE estimator.

``` r
data <- dat2
formula.om <- list()
formula.om[[2]] <- y2 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
  x5 + gam::s(y1)
formula.om[[1]] <- y1 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5

formula.ps <- list()
formula.ps[[2]] <- a ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
  x5  + gam::s(y1)
formula.ps[[1]] <- a ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5

formula.rp <- list()
formula.rp[[2]] <- r2 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
  x5 + gam::s(y1)
formula.rp[[1]] <- r1 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5

formula.pm <- list()
formula.pm[[2]] <- y2 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
  x5 + gam::s(y1)
formula.pm[[1]] <- y1 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5

# create a list of matrices for calibration
k <- 5
inter_mat <- matrix(0, nrow = nrow(data), ncol = choose(k-1, 2)+k-1)
count <- 1
for(i in 1:(k-1)){
  for(j in i:(k-1)){
    inter_mat[,count] <- data[,k+i]*data[,k+j]
    count <- count + 1
  }
}
list.cal <- list()
list.cal[[1]] <- cbind(data.matrix(data[,k + 1:k]), inter_mat) # A
list.cal[[2]] <- cbind(cbind(list.cal[[1]], data$a)) # R1
list.cal[[3]] <- cbind(list.cal[[2]], data$y1,
                       data.matrix(data[,k + 1:(k-1)])*data$y1, data$y1^2) # R2

res.longi <- jtr.longi(formula.ps = formula.ps,
                       formula.rp = formula.rp,
                       formula.om = formula.om,
                       formula.pm = formula.pm,
                       list.cal = list.cal,
                       data = data,
                       type = c("mr", "mr.cal", "rppm"))
res.longi
#> $est
#>        mr    mr.cal      rppm 
#> 1.1897930 0.9718653 1.2532465 
#> 
#> $ve
#>         mr     mr.cal 
#> 0.03732851 0.02691107
```
