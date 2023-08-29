## For the cross-sectional studies ----
rm(list=ls())
library(ggplot2)
library(nleqslv)
library(parallel)

source("sim_1time.R")
source("showres.R")

n <- 500 # sample size
k <- 5 # covariate dimension
alpha <- rep(0.1, k-1)
gamma <- c(1, 1, 1, 1, 0)

numWorkers <- 16
ranseed <- (1:1000)+1234
sim<-length(ranseed)
B <- 500

res <- mclapply(ranseed,main,mc.cores=numWorkers)

save(res, file = "tr_1time.RData")

figure_output(res)

show_res <- assess_function(res)
show_res

xtable::xtable(show_res$table_mat)

## For the longitudinal studies ----
rm(list=ls())
library(mgcv)
library(nleqslv)
library(ggplot2)
library(parallel)
library(xtable)

source("sim_longi.R")
source("showres.R")

n <- 500 # sample size
k <- 5 # covariate dimension
alpha <- rep(0.1, k-1)
gamma1 <- c(rep(20, k-1),0) # obs prob ~ 76%
gamma2 <- c(rep(1, k), 0.1) # obs prob ~ 47%
B <- 500

numWorkers <- 16
ranseed <- (1:1000)+1234
sim<-length(ranseed)

res <- mclapply(ranseed,main,mc.cores=numWorkers)

save(res, file = "mr.RData")

figure_output(res)
show_res <- assess_function(res)
show_res

xtable::xtable(show_res$table_mat_all)
xtable::xtable(show_res$table_mat_ci)


