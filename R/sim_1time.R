#' Generate the data from a simulated longitudinal study
#'
#' @param n the sample size
#' @param k the dimension of baseline covariates
#' @param alpha coefficients of the propensity score
#' @param gamma coefficients of the response probability
#'
#' @return A data frame
#' @export
sim_dat <- function(n, k, alpha, gamma){
  logit_inv <- function(x){
    exp(x) / (1 + exp(x))
  }

  # (1) Generate X
  x_cts <- matrix(rnorm(n*(k-1), 0.25, 1), n, k-1)
  x_bin <- rbinom(n, size = 1, prob = 0.5)
  x <- cbind(x_cts, x_bin)
  colnames(x) <- paste0("x", 1:k)
  # (2) Get Z (nonlinear transformation of X)
  z_cts <- (x_cts^2 + 2*sin(x_cts) - 1.5)/sqrt(2)
  z_bin <- x_bin
  z <- cbind(z_cts, z_bin)
  colnames(z) <- paste0("z", 1:k)
  # (3) Get A|Z
  logit_a <- apply(z_cts, 1, function(x) sum(alpha*x))
  prob_a <- logit_inv(logit_a)
  a <- sapply(prob_a, function(x) rbinom(1,1,prob = x))
  # table(a)
  # (4) Get R|A, Z
  logit_r <- (2*a-1)*apply(z, 1, function(x) sum(gamma*x))/6
  prob_r <- logit_inv(logit_r)
  r <- sapply(prob_r, function(x) rbinom(1,1,prob = x))
  # (5) Get Y|A, Z
  mu_az <- (2 + a)*rowSums(z)/6
  y <- sapply(mu_az, function(x) rnorm(1, x, 1))
  y <- ifelse(y*r==0, NA, y)

  dat_mat <- cbind(x, z, a, r, y)
  dat <- data.frame(dat_mat)

  ry <- ifelse(is.na(r*y), 0, y)
  mu_1z <- (2 + 1)*rowSums(z)/6
  mu_0z <- (2 + 0)*rowSums(z)/6
  logit_pi1 <- (2*1-1)*apply(z, 1, function(x) sum(gamma*x))/6
  logit_pi0 <- (2*0-1)*apply(z, 1, function(x) sum(gamma*x))/6
  pi_1 <- logit_inv(logit_pi1)
  pi_0 <- logit_inv(logit_pi0)
  true_value_psom <- mean(a/prob_a*(ry + (1 - r)*mu_1z) - (1 - a)/(1 - prob_a)*(ry + (1 - r)*mu_0z))
  true_value_rpom <- mean(pi_1*(mu_1z - mu_0z))
  true_value_psrp <- mean(a/prob_a*ry - (1 - a)/(1 - prob_a)*pi_1/pi_0*ry)

  return(list(dat = dat,
              # propensity score output
              prob_a = prob_a,
              # observed probability
              pi_1 = pi_1,
              pi_0 = pi_0,
              # true value
              true_value_psom = true_value_psom,
              true_value_rpom = true_value_rpom,
              true_value_psrp = true_value_psrp))
}

#' Obtain the proposed estimator in the cross-sectional study
#'
#' @param dat the simulated dataset
#' @param k the dimension of baseline covariates
#'
#' @return The point estimates of the eight proposed estimators
#' @import nleqslv
#' @export
model_est <- function(dat, k){
  dat_ctl <- dat
  dat_trt <- dat
  dat_trt$a <- rep(1, nrow(dat))
  dat_ctl$a <- rep(0, nrow(dat))

  # correctly specified (regress on Z)
  ## (c-1) om
  om_fit_c <- lm(y ~ z1 + z2 + z3 + z4 + z5 + factor(a) +
                   factor(a):(z1 + z2 + z3 + z4 + z5),
                 data = dat, na.action = na.omit)
  mu1_c <- predict(om_fit_c, newdata = dat_trt)
  mu0_c <- predict(om_fit_c, newdata = dat_ctl)
  ## (c-2) ps
  ps_fit_c <- glm(a ~ z1 + z2 + z3 + z4 + z5, family = binomial,
                  data = dat)
  e_c <- predict(ps_fit_c, newdata = dat, type = "response")
  ## (c-3) rp
  rp_fit_c <- glm(r ~ z1 + z2 + z3 + z4 + z5 + factor(a) +
                    factor(a):(z1 + z2 + z3 + z4 + z5),
                  family = binomial,
                  data = dat)
  pi1_c <- predict(rp_fit_c, newdata = dat_trt, type = "response")
  pi0_c <- predict(rp_fit_c, newdata = dat_ctl, type = "response")

  # misspecified (regress on X)
  ## (m-1) om
  om_fit_m <- lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(a) +
                   factor(a):(x1 + x2 + x3 + x4 + x5),
                 data = dat, na.action = na.omit)
  mu1_m <- predict(om_fit_m, newdata = dat_trt)
  mu0_m <- predict(om_fit_m, newdata = dat_ctl)
  ## (m-2) ps
  ps_fit_m <- glm(a ~ x1 + x2 + x3 + x4 + x5, family = binomial,
                  data = dat)
  e_m <- predict(ps_fit_m, newdata = dat, type = "response")
  ## (m-3) rp
  rp_fit_m <- glm(r ~ x1 + x2 + x3 + x4 + x5 + factor(a) +
                    factor(a):(x1 + x2 + x3 + x4 + x5),
                  family = binomial,
                  data = dat)
  pi1_m <- predict(rp_fit_m, newdata = dat_trt, type = "response")
  pi0_m <- predict(rp_fit_m, newdata = dat_ctl, type = "response")

  inter_mat <- matrix(0, nrow = n, ncol = k*k)
  for(i in 1:k){
    for(j in 1:k){
      inter_mat[,i + k*(j-1)] <- dat[,k+i]*dat[,k+j]
    }
  }
  g_whole_mat <- data.matrix(dat[,k + 1:k])

  tilde_g <- colMeans(g_whole_mat)

  g1_mat <- g_whole_mat[which(dat$a == 1),]
  g0_mat <- g_whole_mat[which(dat$a == 0),]

  weight_fn <- function(g_mat, tilde_g){
    opt_fn <- function(lambda){
      first_part <- as.vector(exp(g_mat%*%lambda) + 1)/sum(exp(g_mat%*%lambda) + 1)
      colSums(first_part*g_mat) - tilde_g
    }
    set.seed(123)
    ini_lambda <- matrix(rnorm(length(tilde_g)*10,0,0.5),nrow = 10)
    temp <- nleqslv::searchZeros(ini_lambda, opt_fn)
    ind_temp <- 1

    w <- temp$x[ind_temp,]
    deno <- apply(g_mat, 1, function(x) exp(sum(w*x)))
    num <- sum(deno)
    res <- deno/num
    return(res)
  }

  w1_cal <- weight_fn(g1_mat, tilde_g)
  w0_cal <- weight_fn(g0_mat, tilde_g)

  gr_mat <- g_whole_mat[which(dat$r == 1),]
  wr_cal <- weight_fn(gr_mat, tilde_g)
  ind_temp1 <- which(dat[which(dat$a == 0), ]$r == 1)
  ind_temp2 <- which(dat[which(dat$r == 1), ]$a == 0)
  w0r_cal <- (w0_cal[ind_temp1]*wr_cal[ind_temp2])/sum(w0_cal[ind_temp1]*wr_cal[ind_temp2])

  rpom_est <- function(pi1, mu1, mu0){
    return(mean(pi1*(mu1 - mu0)))
  }

  psom_est <- function(e, mu1, mu0){
    ry <- dat$r*dat$y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part <- dat$a/e*(ry + (1 - dat$r)*mu0)
    second_part <- (1 - dat$a)/(1 - e)*(ry + (1 - dat$r)*mu0)
    return(mean(first_part - second_part))
  }

  psom_alter_est <- function(e, mu1, mu0){
    ry <- dat$r*dat$y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part_up <- dat$a/e*(ry + (1 - dat$r)*mu0)
    first_part_down <- dat$a/e
    second_part_up <- (1 - dat$a)/(1 - e)*(ry + (1 - dat$r)*mu0)
    second_part_down <- (1 - dat$a)/(1 - e)
    return(sum(first_part_up)/sum(first_part_down) - sum(second_part_up)/sum(second_part_down))
  }

  psrp_est <- function(e, pi1, pi0){
    ry <- dat$r*dat$y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part <- dat$a/e*ry
    second_part <- (1 - dat$a)/(1 - e)*(pi1/pi0*ry)
    return(mean(first_part - second_part))
  }

  psrp_alter_est <- function(e, pi1, pi0){
    ry <- dat$r*dat$y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part_up <- dat$a/e*ry
    first_part_down <- dat$a/e
    second_part_up <- (1 - dat$a)/(1 - e)*(pi1/pi0*ry)
    second_part_down <- (1 - dat$a)*dat$r/((1 - e)*pi0)
    return(sum(first_part_up)/sum(first_part_down) - sum(second_part_up)/sum(second_part_down))
  }

  psrpom_est <- function(e, pi1, pi0, mu1, mu0){
    ry_mu0 <- dat$r*(dat$y - mu0)
    ry_mu0 <- ifelse(is.na(ry_mu0), 0, ry_mu0)
    first_part <- (dat$a/e - (1 - dat$a)/(1 - e)*pi1/pi0)*ry_mu0
    second_part <- (dat$a - e)/e*pi1*(mu1 - mu0)
    point_value <- mean(first_part - second_part)
    var_value <- mean((first_part - second_part - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }

  psrpom_alter_est <- function(e, pi1, pi0, mu1, mu0){
    ry_mu0 <- dat$r*(dat$y - mu0)
    ry_mu0 <- ifelse(is.na(ry_mu0), 0, ry_mu0)
    first_part_up <- dat$a/e*(ry_mu0 - pi1*(mu1 - mu0))
    first_part_down <- dat$a/e
    second_part_up <- (1 - dat$a)/(1 - e)*pi1/pi0*ry_mu0
    second_part_down <- (1 - dat$a)*dat$r/((1 - e)*pi0)
    third_part <- pi1*(mu1 - mu0)
    point_value <- sum(first_part_up)/sum(first_part_down) -
      sum(second_part_up)/sum(second_part_down) + mean(third_part)
    var_value <- mean((first_part_up/mean(first_part_down) -
      second_part_up/mean(second_part_down) + third_part - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }

  psrpom_weight_est <- function(w1, w0r, pi1, mu1, mu0){
    ry_mu0 <- dat$r*(dat$y - mu0)
    ry_mu0 <- ifelse(is.na(ry_mu0), 0, ry_mu0)
    diff <- ry_mu0 - pi1*(mu1 - mu0)
    diff1 <- diff[which(dat$a == 1)]
    piry_mu00r <- (pi1*ry_mu0)[which(dat$a == 0 & dat$r == 1)]
    first_part <- sum(w1*diff1)
    second_part <- sum(w0r*piry_mu00r)
    third_part <- mean(pi1*(mu1 - mu0))
    point_value <- first_part - second_part + third_part

    part1_long <- rep(0, n)
    part1_long[dat$a == 1]<- n*w1*diff1
    part2_long <- rep(0, n)
    part2_long[dat$a == 0 & dat$r == 1]<- n*w0r*piry_mu00r
    part3_long <- pi1*(mu1 - mu0)

    var_value <- mean((part1_long - part2_long + part3_long - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }

  # CASE 1: c(ps, om, rp)
  tr_res <- psrpom_est(e_c, pi1_c, pi0_c, mu1_c, mu0_c)
  tr_psrpom_psrpom_est <- tr_res$point_value
  tr_psrpom_psrpom_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_c, pi1_c, pi0_c, mu1_c, mu0_c)
  tr_psrpom_psrpom_alter_est <- tr_alter_res$point_value
  tr_psrpom_psrpom_alter_var <- tr_alter_res$var_value
  psrp_psrpom_est <- psrp_est(e_c, pi1_c, pi0_c)
  psrp_psrpom_alter_est <- psrp_alter_est(e_c, pi1_c, pi0_c)
  psom_psrpom_est <- psom_est(e_c, mu1_c, mu0_c)
  psom_psrpom_alter_est <- psom_alter_est(e_c, mu1_c, mu0_c)
  rpom_psrpom_est <- rpom_est(pi1_c, mu1_c, mu0_c)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_c, mu1_c, mu0_c)
  tr_cal_psrpom_est <- tr_cal_res$point_value
  tr_cal_psrpom_var <- tr_cal_res$var_value

  psrpom_all <- c(tr_psrpom_psrpom_est,
                  tr_psrpom_psrpom_alter_est,
                  tr_cal_psrpom_est,
                  psrp_psrpom_est,
                  psrp_psrpom_alter_est,
                  psom_psrpom_est,
                  psom_psrpom_alter_est,
                  rpom_psrpom_est)
  psrpom_var <- c(tr_psrpom_psrpom_var,
                  tr_psrpom_psrpom_alter_var,
                  tr_cal_psrpom_var)

  # CASE 2: c(ps, om); m(rp)
  tr_res <- psrpom_est(e_c, pi1_m, pi0_m, mu1_c, mu0_c)
  tr_psrpom_psom_est <- tr_res$point_value
  tr_psrpom_psom_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_c, pi1_m, pi0_m, mu1_c, mu0_c)
  tr_psrpom_psom_alter_est <- tr_alter_res$point_value
  tr_psrpom_psom_alter_var <- tr_alter_res$var_value
  psrp_psom_est <- psrp_est(e_c, pi1_m, pi0_m)
  psrp_psom_alter_est <- psrp_alter_est(e_c, pi1_m, pi0_m)
  psom_psom_est <- psom_est(e_c, mu1_c, mu0_c)
  psom_psom_alter_est <- psom_alter_est(e_c, mu1_c, mu0_c)
  rpom_psom_est <- rpom_est(pi1_m, mu1_c, mu0_c)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_m, mu1_c, mu0_c)
  tr_cal_psom_est <- tr_cal_res$point_value
  tr_cal_psom_var <- tr_cal_res$var_value

  psom_all <- c(tr_psrpom_psom_est,
                tr_psrpom_psom_alter_est,
                tr_cal_psom_est,
                psrp_psom_est,
                psrp_psom_alter_est,
                psom_psom_est,
                psom_psom_alter_est,
                rpom_psom_est)
  psom_var <- c(tr_psrpom_psom_var,
                tr_psrpom_psom_alter_var,
                tr_cal_psom_var)

  # CASE 3: c(ps); m(om, rp)
  tr_res <- psrpom_est(e_c, pi1_m, pi0_m, mu1_m, mu0_m)
  tr_psrpom_ps_est <- tr_res$point_value
  tr_psrpom_ps_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_c, pi1_m, pi0_m, mu1_m, mu0_m)
  tr_psrpom_ps_alter_est <- tr_alter_res$point_value
  tr_psrpom_ps_alter_var <- tr_alter_res$var_value
  psrp_ps_est <- psrp_est(e_c, pi1_m, pi0_m)
  psrp_ps_alter_est <- psrp_alter_est(e_c, pi1_m, pi0_m)
  psom_ps_est <- psom_est(e_c, mu1_m, mu0_m)
  psom_ps_alter_est <- psom_alter_est(e_c, mu1_m, mu0_m)
  rpom_ps_est <- rpom_est(pi1_m, mu1_m, mu0_m)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_m, mu1_m, mu0_m)
  tr_cal_ps_est <- tr_cal_res$point_value
  tr_cal_ps_var <- tr_cal_res$var_value

  ps_all <- c(tr_psrpom_ps_est,
              tr_psrpom_ps_alter_est,
              tr_cal_ps_est,
              psrp_ps_est,
              psrp_ps_alter_est,
              psom_ps_est,
              psom_ps_alter_est,
              rpom_ps_est)
  ps_var <- c(tr_psrpom_ps_var,
              tr_psrpom_ps_alter_var,
              tr_cal_ps_var)

  # CASE 4: c(rp, om); m(ps)
  tr_res <- psrpom_est(e_m, pi1_c, pi0_c, mu1_c, mu0_c)
  tr_psrpom_rpom_est <- tr_res$point_value
  tr_psrpom_rpom_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_m, pi1_c, pi0_c, mu1_c, mu0_c)
  tr_psrpom_rpom_alter_est <- tr_alter_res$point_value
  tr_psrpom_rpom_alter_var <- tr_alter_res$var_value
  psrp_rpom_est <- psrp_est(e_m, pi1_c, pi0_c)
  psrp_rpom_alter_est <- psrp_alter_est(e_m, pi1_c, pi0_c)
  psom_rpom_est <- psom_est(e_m, mu1_c, mu0_c)
  psom_rpom_alter_est <- psom_alter_est(e_m, mu1_c, mu0_c)
  rpom_rpom_est <- rpom_est(pi1_c, mu1_c, mu0_c)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_c, mu1_c, mu0_c)
  tr_cal_rpom_est <- tr_cal_res$point_value
  tr_cal_rpom_var <- tr_cal_res$var_value

  rpom_all <- c(tr_psrpom_rpom_est,
                tr_psrpom_rpom_alter_est,
                tr_cal_rpom_est,
                psrp_rpom_est,
                psrp_rpom_alter_est,
                psom_rpom_est,
                psom_rpom_alter_est,
                rpom_rpom_est)
  rpom_var <- c(tr_psrpom_rpom_var,
                tr_psrpom_rpom_alter_var,
                tr_cal_rpom_var)

  # CASE 5: c(rp); m(ps, om)
  tr_res <- psrpom_est(e_m, pi1_c, pi0_c, mu1_m, mu0_m)
  tr_psrpom_rp_est <- tr_res$point_value
  tr_psrpom_rp_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_m, pi1_c, pi0_c, mu1_m, mu0_m)
  tr_psrpom_rp_alter_est <- tr_alter_res$point_value
  tr_psrpom_rp_alter_var <- tr_alter_res$var_value
  psrp_rp_est <- psrp_est(e_m, pi1_c, pi0_c)
  psrp_rp_alter_est <- psrp_alter_est(e_m, pi1_c, pi0_c)
  psom_rp_est <- psom_est(e_m, mu1_m, mu0_m)
  psom_rp_alter_est <- psom_alter_est(e_m, mu1_m, mu0_m)
  rpom_rp_est <- rpom_est(pi1_c, mu1_m, mu0_m)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_c, mu1_m, mu0_m)
  tr_cal_rp_est <- tr_cal_res$point_value
  tr_cal_rp_var <- tr_cal_res$var_value
  rp_all <- c(tr_psrpom_rp_est,
              tr_psrpom_rp_alter_est,
              tr_cal_rp_est,
              psrp_rp_est,
              psrp_rp_alter_est,
              psom_rp_est,
              psom_rp_alter_est,
              rpom_rp_est)
  rp_var <- c(tr_psrpom_rp_var,
              tr_psrpom_rp_alter_var,
              tr_cal_rp_var)

  # CASE 6: c(rp, ps); m(om)
  tr_res <- psrpom_est(e_c, pi1_c, pi0_c, mu1_m, mu0_m)
  tr_psrpom_psrp_est <- tr_res$point_value
  tr_psrpom_psrp_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_c, pi1_c, pi0_c, mu1_m, mu0_m)
  tr_psrpom_psrp_alter_est <- tr_alter_res$point_value
  tr_psrpom_psrp_alter_var <- tr_alter_res$var_value
  psrp_psrp_est <- psrp_est(e_c, pi1_c, pi0_c)
  psrp_psrp_alter_est <- psrp_alter_est(e_c, pi1_c, pi0_c)
  psom_psrp_est <- psom_est(e_c, mu1_m, mu0_m)
  psom_psrp_alter_est <- psom_alter_est(e_c, mu1_m, mu0_m)
  rpom_psrp_est <- rpom_est(pi1_c, mu1_m, mu0_m)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_c, mu1_m, mu0_m)
  tr_cal_psrp_est <- tr_cal_res$point_value
  tr_cal_psrp_var <- tr_cal_res$var_value
  psrp_all <- c(tr_psrpom_psrp_est,
                tr_psrpom_psrp_alter_est,
                tr_cal_psrp_est,
                psrp_psrp_est,
                psrp_psrp_alter_est,
                psom_psrp_est,
                psom_psrp_alter_est,
                rpom_psrp_est)
  psrp_var <- c(tr_psrpom_psrp_var,
                tr_psrpom_psrp_alter_var,
                tr_cal_psrp_var)

  # CASE 7: c(om); m(rp, ps)
  tr_res <- psrpom_est(e_m, pi1_m, pi0_m, mu1_c, mu0_c)
  tr_psrpom_om_est <- tr_res$point_value
  tr_psrpom_om_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_m, pi1_m, pi0_m, mu1_c, mu0_c)
  tr_psrpom_om_alter_est <- tr_alter_res$point_value
  tr_psrpom_om_alter_var <- tr_alter_res$var_value
  psrp_om_est <- psrp_est(e_m, pi1_m, pi0_m)
  psrp_om_alter_est <- psrp_alter_est(e_m, pi1_m, pi0_m)
  psom_om_est <- psom_est(e_m, mu1_c, mu0_c)
  psom_om_alter_est <- psom_alter_est(e_m, mu1_c, mu0_c)
  rpom_om_est <- rpom_est(pi1_m, mu1_c, mu0_c)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_m, mu1_c, mu0_c)
  tr_cal_om_est <- tr_cal_res$point_value
  tr_cal_om_var <- tr_cal_res$var_value
  om_all <- c(tr_psrpom_om_est,
              tr_psrpom_om_alter_est,
              tr_cal_om_est,
              psrp_om_est,
              psrp_om_alter_est,
              psom_om_est,
              psom_om_alter_est,
              rpom_om_est)
  om_var <- c(tr_psrpom_om_var,
              tr_psrpom_om_alter_var,
              tr_cal_om_var)

  # CASE 8: m(ps, om, rp)
  tr_res <- psrpom_est(e_m, pi1_m, pi0_m, mu1_m, mu0_m)
  tr_psrpom_none_est <- tr_res$point_value
  tr_psrpom_none_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_m, pi1_m, pi0_m, mu1_m, mu0_m)
  tr_psrpom_none_alter_est <- tr_alter_res$point_value
  tr_psrpom_none_alter_var <- tr_alter_res$var_value
  psrp_none_est <- psrp_est(e_m, pi1_m, pi0_m)
  psrp_none_alter_est <- psrp_alter_est(e_m, pi1_m, pi0_m)
  psom_none_est <- psom_est(e_m, mu1_m, mu0_m)
  psom_none_alter_est <- psom_alter_est(e_m, mu1_m, mu0_m)
  rpom_none_est <- rpom_est(pi1_m, mu1_m, mu0_m)
  tr_cal_res <- psrpom_weight_est(w1_cal, w0r_cal, pi1_m, mu1_m, mu0_m)
  tr_cal_none_est <- tr_cal_res$point_value
  tr_cal_none_var <- tr_cal_res$var_value

  none_all <- c(tr_psrpom_none_est,
                tr_psrpom_none_alter_est,
                tr_cal_none_est,
                psrp_none_est,
                psrp_none_alter_est,
                psom_none_est,
                psom_none_alter_est,
                rpom_none_est)
  none_var <- c(tr_psrpom_none_var,
                tr_psrpom_none_alter_var,
                tr_cal_none_var)

  var_all <- c(psrpom_var,
               psom_var,
               ps_var,
               rpom_var,
               rp_var,
               psrp_var,
               om_var,
               none_var)

  return(list(psrpom_all = psrpom_all,
              psom_all = psom_all,
              ps_all = ps_all,
              rpom_all = rpom_all,
              rp_all = rp_all,
              psrp_all = psrp_all,
              om_all = om_all,
              none_all = none_all,
              psrpom_var = psrpom_var,
              psom_var = psom_var,
              ps_var = ps_var,
              rpom_var = rpom_var,
              rp_var = rp_var,
              psrp_var = psrp_var,
              om_var = om_var,
              none_var = none_var))
}

#' Nonparametric bootstrap to obtain the variance estimates
#'
#' @param dat the simulated dataset
#' @param B the number of bootstrap replicates
#' @param psrpom_est a vector of the point estimates of the eight proposed
#' estimators when three models are correctly specified
#' @param psom_est a vector of the point estimates of the eight proposed
#' estimators when (ps, om) are correctly specified
#' @param ps_est a vector of the point estimates of the eight proposed
#' estimators when (ps) is correctly specified
#' @param rpom_est a vector of the point estimates of the eight proposed
#' estimators when (rp, om) are correctly specified
#' @param rp_est a vector of the point estimates of the eight proposed
#' estimators when (rp) is correctly specified
#' @param psrp_est a vector of the point estimates of the eight proposed
#' estimators when (ps, rp) are correctly specified
#' @param om_est a vector of the point estimates of the eight proposed
#' estimators when (om) is correctly specified
#' @param none_est a vector of the point estimates of the eight proposed
#' estimators when none of the modesl are correctly specified
#'
#' @return the estimated variation of the estimators using nonparmetric
#' bootstrap, symmetric-t bootstrap, and bootstrap percentile
#' @export
nonpara_fn <- function(dat, B, psrpom_est, psom_est, ps_est,
                       rpom_est, rp_est, psrp_est, om_est, none_est){
  n <- nrow(dat)
  psrpom_all_np <- matrix(0, 8, B)
  psom_all_np <- matrix(0, 8, B)
  ps_all_np <- matrix(0, 8, B)
  rpom_all_np <- matrix(0, 8, B)
  rp_all_np <- matrix(0, 8, B)
  psrp_all_np <- matrix(0, 8, B)
  om_all_np <- matrix(0, 8, B)
  none_all_np <- matrix(0, 8, B)
  psrpom_var <- matrix(0, 3, B)
  psom_var <- matrix(0, 3, B)
  ps_var <- matrix(0, 3, B)
  rpom_var <- matrix(0, 3, B)
  rp_var <- matrix(0, 3, B)
  psrp_var <- matrix(0, 3, B)
  om_var <- matrix(0, 3, B)
  none_var <- matrix(0, 3, B)
  t_psrpom_star <- matrix(0, 3, B)
  t_psom_star <- matrix(0, 3, B)
  t_ps_star <- matrix(0, 3, B)
  t_rpom_star <- matrix(0, 3, B)
  t_rp_star <- matrix(0, 3, B)
  t_psrp_star <- matrix(0, 3, B)
  t_om_star <- matrix(0, 3, B)
  t_none_star <- matrix(0, 3, B)

  for(b in 1:B){
    set.seed(b)
    boot_id <- sample(1:n, size = n, replace = TRUE)
    boot_dat <- dat[boot_id,]
    boot_res <- model_est(boot_dat)
    psrpom_all_np[,b] <- boot_res$psrpom_all
    psom_all_np[,b] <- boot_res$psom_all
    ps_all_np[,b] <- boot_res$ps_all
    rpom_all_np[,b] <- boot_res$rpom_all
    rp_all_np[,b] <- boot_res$rp_all
    psrp_all_np[,b] <- boot_res$psrp_all
    om_all_np[,b] <- boot_res$om_all
    none_all_np[,b] <- boot_res$none_all

    psrpom_var[,b] <- boot_res$psrpom_var
    psom_var[,b] <- boot_res$psom_var
    ps_var[,b] <- boot_res$ps_var
    rpom_var[,b] <- boot_res$rpom_var
    rp_var[,b] <- boot_res$rp_var
    psrp_var[,b] <- boot_res$psrp_var
    om_var[,b] <- boot_res$om_var
    none_var[,b] <- boot_res$none_var

    # compute t-star in each case (total: 8)
    t_psrpom_star[,b] <- (psrpom_all_np[1:3,b] - psrpom_est)/sqrt(psrpom_var[1:3,b])
    t_psom_star[,b] <- (psom_all_np[1:3,b] - psom_est)/sqrt(psom_var[1:3,b])
    t_ps_star[,b] <- (ps_all_np[1:3,b] - ps_est)/sqrt(ps_var[1:3,b])
    t_rpom_star[,b] <- (rpom_all_np[1:3,b] - rpom_est)/sqrt(rpom_var[1:3,b])
    t_rp_star[,b] <- (rp_all_np[1:3,b] - rp_est)/sqrt(rp_var[1:3,b])
    t_psrp_star[,b] <- (psrp_all_np[1:3,b] - psrp_est)/sqrt(psrp_var[1:3,b])
    t_om_star[,b] <- (om_all_np[1:3,b] - om_est)/sqrt(om_var[1:3,b])
    t_none_star[,b] <- (none_all_np[1:3,b] - none_est)/sqrt(none_var[1:3,b])
  }
  var_psrpom_all_np <- apply(psrpom_all_np, 1, var)
  var_psom_all_np <- apply(psom_all_np, 1, var)
  var_ps_all_np <- apply(ps_all_np, 1, var)
  var_rpom_all_np <- apply(rpom_all_np, 1, var)
  var_rp_all_np <- apply(rp_all_np, 1, var)
  var_psrp_all_np <- apply(psrp_all_np, 1, var)
  var_om_all_np <- apply(om_all_np, 1, var)
  var_none_all_np <- apply(none_all_np, 1, var)

  # compute c-star in each case
  c_psrpom_star <- apply(abs(t_psrpom_star), 1, function(x) quantile(x, 0.95))
  c_psom_star <- apply(abs(t_psom_star), 1, function(x) quantile(x, 0.95))
  c_ps_star <- apply(abs(t_ps_star), 1, function(x) quantile(x, 0.95))
  c_rpom_star <- apply(abs(t_rpom_star), 1, function(x) quantile(x, 0.95))
  c_rp_star <- apply(abs(t_rp_star), 1, function(x) quantile(x, 0.95))
  c_psrp_star <- apply(abs(t_psrp_star), 1, function(x) quantile(x, 0.95))
  c_om_star <- apply(abs(t_om_star), 1, function(x) quantile(x, 0.95))
  c_none_star <- apply(abs(t_none_star), 1, function(x) quantile(x, 0.95))

  # percentile bootstrap
  percentile_psrpom_boot <- apply(psrpom_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_psom_boot <- apply(psom_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_ps_boot <- apply(ps_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_rpom_boot <- apply(rpom_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_rp_boot <- apply(rp_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_psrp_boot <- apply(psrp_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_om_boot <- apply(om_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  percentile_none_boot <- apply(none_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))

  return(list(var_psrpom_all_np = var_psrpom_all_np,
              var_psom_all_np = var_psom_all_np,
              var_ps_all_np = var_ps_all_np,
              var_rpom_all_np = var_rpom_all_np,
              var_rp_all_np = var_rp_all_np,
              var_psrp_all_np = var_psrp_all_np,
              var_om_all_np = var_om_all_np,
              var_none_all_np = var_none_all_np,
              # c-star for symmetric t bootstrap
              c_psrpom_star = c_psrpom_star,
              c_psom_star = c_psom_star,
              c_ps_star = c_ps_star,
              c_rpom_star = c_rpom_star,
              c_rp_star = c_rp_star,
              c_psrp_star = c_psrp_star,
              c_om_star = c_om_star,
              c_none_star = c_none_star,
              # percentile bootstrap
              percentile_psrpom_boot = percentile_psrpom_boot,
              percentile_psom_boot = percentile_psom_boot,
              percentile_ps_boot = percentile_ps_boot,
              percentile_rpom_boot = percentile_rpom_boot,
              percentile_rp_boot = percentile_rp_boot,
              percentile_psrp_boot = percentile_psrp_boot,
              percentile_om_boot = percentile_om_boot,
              percentile_none_boot = percentile_none_boot
              ))
}


#' Main function to get both point and variance estimates
#'
#' @param seed a specific seed number
#' @param n the sample size
#' @param k the dimension of baseline covariates
#' @param alpha coefficients of the propensity score
#' @param gamma coefficients of the response probability
#' @param B the number of bootstrap replicates
#'
#' @return the estimated point and variance estimates
#' @export
main <- function(seed, n, k, alpha, gamma, B){
  set.seed(seed)
  dat_list <- sim_dat(n = n, k = k, alpha = alpha, gamma = gamma)
  dat <- dat_list$dat
  true_value_psom <- dat_list$true_value_psom
  true_value_rpom <- dat_list$true_value_rpom
  true_value_psrp <- dat_list$true_value_psrp
  prob_a <- dat_list$prob_a # monitor the propensity score
  pi_1 <- dat_list$pi_1
  pi_0 <- dat_list$pi_0
  est_value <- model_est(dat, k)
  psrpom_all <- est_value$psrpom_all
  psom_all <- est_value$psom_all
  ps_all <- est_value$ps_all
  rpom_all <- est_value$rpom_all
  rp_all <- est_value$rp_all
  psrp_all <- est_value$psrp_all
  om_all <- est_value$om_all
  none_all <- est_value$none_all
  w1_cal <- est_value$w1_cal
  w0_cal <- est_value$w0_cal
  # variance estimation (asymptotic theory)
  var_psrpom_all_asym <- est_value$psrpom_var
  var_psom_all_asym <- est_value$psom_var
  var_ps_all_asym <- est_value$ps_var
  var_rpom_all_asym <- est_value$rpom_var
  var_rp_all_asym <- est_value$rp_var
  var_psrp_all_asym <- est_value$psrp_var
  var_om_all_asym <- est_value$om_var
  var_none_all_asym <- est_value$none_var

  # variance estimation (bootstrap)
  var_np <- nonpara_fn(dat, B, psrpom_est = psrpom_all[1:3],
                       psom_est = psom_all[1:3],
                       ps_est = ps_all[1:3],
                       rpom_est = rpom_all[1:3],
                       rp_est = rp_all[1:3],
                       psrp_est = psrp_all[1:3],
                       om_est = om_all[1:3],
                       none_est = none_all[1:3])
  var_psrpom_all_np <- var_np$var_psrpom_all_np
  var_psom_all_np <- var_np$var_psom_all_np
  var_ps_all_np <- var_np$var_ps_all_np
  var_rpom_all_np <- var_np$var_rpom_all_np
  var_rp_all_np <- var_np$var_rp_all_np
  var_psrp_all_np <- var_np$var_psrp_all_np
  var_om_all_np <- var_np$var_om_all_np
  var_none_all_np <- var_np$var_none_all_np

  c_psrpom_star <- var_np$c_psrpom_star
  c_psom_star <- var_np$c_psom_star
  c_ps_star <- var_np$c_ps_star
  c_rpom_star <- var_np$c_rpom_star
  c_rp_star <- var_np$c_rp_star
  c_psrp_star <- var_np$c_psrp_star
  c_om_star <- var_np$c_om_star
  c_none_star <- var_np$c_none_star

  percentile_psrpom_boot <- var_np$percentile_psrpom_boot
  percentile_psom_boot <- var_np$percentile_psom_boot
  percentile_ps_boot <- var_np$percentile_ps_boot
  percentile_rpom_boot <- var_np$percentile_rpom_boot
  percentile_rp_boot <- var_np$percentile_rp_boot
  percentile_psrp_boot <- var_np$percentile_psrp_boot
  percentile_om_boot <- var_np$percentile_om_boot
  percentile_none_boot <- var_np$percentile_none_boot

  return(list(psrpom_all = psrpom_all,
              psom_all = psom_all,
              ps_all = ps_all,
              rpom_all = rpom_all,
              rp_all = rp_all,
              psrp_all = psrp_all,
              om_all = om_all,
              none_all = none_all,
              # true value
              true_value_psom = true_value_psom,
              true_value_rpom = true_value_rpom,
              true_value_psrp = true_value_psrp,
              # variance estimate (nonparametric boot)
              var_psrpom_all_np = var_psrpom_all_np,
              var_psom_all_np = var_psom_all_np,
              var_ps_all_np = var_ps_all_np,
              var_rpom_all_np = var_rpom_all_np,
              var_rp_all_np = var_rp_all_np,
              var_psrp_all_np = var_psrp_all_np,
              var_om_all_np = var_om_all_np,
              var_none_all_np = var_none_all_np,
              # variance estimate (asymptotic theory)
              var_psrpom_all_asym = var_psrpom_all_asym,
              var_psom_all_asym = var_psom_all_asym,
              var_ps_all_asym = var_ps_all_asym,
              var_rpom_all_asym = var_rpom_all_asym,
              var_rp_all_asym = var_rp_all_asym,
              var_psrp_all_asym = var_psrp_all_asym,
              var_om_all_asym = var_om_all_asym,
              var_none_all_asym = var_none_all_asym,
              # c-star
              c_psrpom_star = c_psrpom_star,
              c_psom_star = c_psom_star,
              c_ps_star = c_ps_star,
              c_rpom_star = c_rpom_star,
              c_rp_star = c_rp_star,
              c_psrp_star = c_psrp_star,
              c_om_star = c_om_star,
              c_none_star = c_none_star,
              # percentile bootstrap (CI)
              percentile_psrpom_boot = percentile_psrpom_boot,
              percentile_psom_boot = percentile_psom_boot,
              percentile_ps_boot = percentile_ps_boot,
              percentile_rpom_boot = percentile_rpom_boot,
              percentile_rp_boot = percentile_rp_boot,
              percentile_psrp_boot = percentile_psrp_boot,
              percentile_om_boot = percentile_om_boot,
              percentile_none_boot = percentile_none_boot
              ))
}

#' @importFrom stats binomial gaussian na.omit predict quantile rbinom rnorm var
NULL
