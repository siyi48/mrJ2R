sim_dat <- function(n, k, alpha, gamma1, gamma2){
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
  # (4) Get R1|A, Z
  logit_r1 <- (2*a-1)*apply(z, 1, function(x) sum(gamma1*x))/36
  prob_r1 <- logit_inv(logit_r1)
  r1 <- sapply(prob_r1, function(x) rbinom(1,1,prob = x))
  # (5) Get Y1|A, Z
  mu_az <- (2 + a)*rowSums(z)/6
  y1 <- sapply(mu_az, function(x) rnorm(1, x, 1))
  y1 <- ifelse(y1*r1==0, NA, y1)
  # (6) Get R2|A, Z, Y1 (need to reset parameters)
  # only for observed data
  r2 <- r1
  logit_r2 <- (2*a[r1==1]-1)*apply(cbind(z, y1)[r1==1,], 1, function(x) sum(gamma2*x))/6
  prob_r2 <- logit_inv(logit_r2)
  r2[r1==1] <- sapply(prob_r2, function(x) rbinom(1,1,prob = x))
  # (7) Get Y2|A, Z, Y1
  y2 <- y1
  mu_azy1 <- (2 + a[r1==1])*rowSums(cbind(z, y1)[r1==1,])/3
  y2[r1 == 1] <- sapply(mu_azy1, function(x) rnorm(1, x, 1))
  y2 <- ifelse(y2*r2==0, NA, y2)
  
  dat_mat <- cbind(x, z, a, r1, r2, y1, y2)
  dat <- data.frame(dat_mat)
  
  # true model
  y2_imp <- y2
  mu_0z <- 2*rowSums(z)/6
  mu_0zy1 <- 2*rowSums(cbind(z, y1)[r1==1,])/3
  y2_imp[r1 == 1 & r2 == 0] <- mu_0zy1[which(r2[r1==1] == 0)]
  y2_imp[r1 == 0] <- 2*rowSums(cbind(z, mu_0z)[r1 == 0,])/3
  long_mu_0zy1 <- r2
  long_mu_0zy1[r1 == 1] <- mu_0zy1
  r2y2 <- ifelse(is.na(r2*y2), 0, y2)
  true_value_psom <- mean(a/prob_a*(r2y2 + r1*(1 - r2)*long_mu_0zy1 + 
                                      (1 - r1)*2*rowSums(cbind(z, mu_0z))/3) - 
                            (1 - a)/(1 - prob_a)*(r2y2 + r1*(1 - r2)*long_mu_0zy1 + 
                                                    (1 - r1)*2*rowSums(cbind(z, mu_0z))/3))
 
  return(list(dat = dat, 
              # propensity score output
              prob_a = prob_a,
              # observed probability
              # pi_1 = pi_1,
              # pi_0 = pi_0,
              # # true value
              true_value_psom = true_value_psom
  ))
}

model_est <- function(dat){
  dat_ctl <- dat
  dat_trt <- dat
  dat_trt$a <- rep(1, nrow(dat))
  dat_ctl$a <- rep(0, nrow(dat))
  dat_a1 <- dat[dat$a == 1,]
  dat_a0 <- dat[dat$a == 0,]
  dat_R11 <- dat[dat$r1 == 1,]
  dat_R11a1 <- dat_R11[dat_R11$a == 1,]
  dat_R11a0 <- dat_R11[dat_R11$a == 0,]
  dat_ctlR11 <- dat_ctl[dat_ctl$r1 == 1,]
  dat_trtR11 <- dat_trt[dat_trt$r1 == 1,]
  
  ## (c-1) om
  # second time point
  om2_fit_a1 <- gam(y2 ~ s(x1) + s(x2) + s(x3) + 
                      s(x4) + x5 + s(y1),
                    data = dat_R11a1, na.action = na.omit)
  om2_fit_a0 <- gam(y2 ~ s(x1) + s(x2) + s(x3) + s(x4) + x5 + s(y1),
                    data = dat_R11a0, na.action = na.omit)
  mu21_c <- predict(om2_fit_a1, newdata = dat_R11)
  mu20_c <- predict(om2_fit_a0, newdata = dat_R11)
  
  # first time point (backward)
  om1_fit_a0 <- gam(mu20_c[dat_R11$a == 0] ~ s(x1) + s(x2) + s(x3) + s(x4) + x5,
                    data = dat_R11a0, na.action = na.omit)
  mu10_c <- predict(om1_fit_a0, newdata = dat_ctl)
  
  ## (c-2) ps (keep the same)
  # ps_fit_c <- gam(a ~ s(z1) + s(z2) + s(z3) + s(z4) + z5, 
  #                 family = binomial,
  #                 data = dat)
  ps_fit_c <- gam(a ~ s(x1) + s(x2) + s(x3) + s(x4) + x5, 
                  family = binomial,
                  data = dat)
  e_c <- predict(ps_fit_c, newdata = dat, type = "response")
  
  # for a weird ratio
  # ps2_fit <- gam(a ~ s(z1) + s(z2) + s(z3) + s(z4) + z5 + s(y1),
  #                family = binomial,
  #                data = dat_R11)
  ps2_fit <- gam(a ~ s(x1) + s(x2) + s(x3) + s(x4) + x5 + s(y1),
                 family = binomial,
                 data = dat_R11)
  eh1_c <- predict(ps2_fit, newdata = dat_R11, type = "response")
  # ps3_fit <- gam(a ~ s(z1) + s(z2) + s(z3) + s(z4) + z5,
  #                family = binomial,
  #                data = dat_R11)
  ps3_fit <- gam(a ~ s(x1) + s(x2) + s(x3) + s(x4) + x5,
                 family = binomial,
                 data = dat_R11)
  eh0_c <- predict(ps3_fit, newdata = dat_R11, type = "response")
  ps_ratio_c <- (eh1_c/eh0_c)/((1 - eh1_c)/(1 - eh0_c))

  ## (c-3) rp
  # second time point
  # rp2_fit_a1 <- gam(r2 ~ s(z1) + s(z2) + s(z3) + s(z4) + z5 + s(y1),
  #                   family = binomial,
  #                   data = dat_R11a1)
  # rp2_fit_a0 <- gam(r2 ~ s(z1) + s(z2) + s(z3) + s(z4) + z5 + s(y1),
  #                   family = binomial,
  #                   data = dat_R11a0)
  rp2_fit_a1 <- gam(r2 ~ s(x1) + s(x2) + s(x3) + s(x4) + x5 + s(y1),
                    family = binomial,
                    data = dat_R11a1)
  rp2_fit_a0 <- gam(r2 ~ s(x1) + s(x2) + s(x3) + s(x4) + x5 + s(y1),
                    family = binomial,
                    data = dat_R11a0)
  pi21_c <- predict(rp2_fit_a1, newdata = dat_trtR11, type = "response")
  pi20_c <- predict(rp2_fit_a0, newdata = dat_ctlR11, type = "response")
  
  # first time point
  # rp1_fit_a1 <- gam(r1 ~ s(z1) + s(z2) + s(z3) + s(z4) + z5,
  #                   family = binomial,
  #                   data = dat_a1)
  # rp1_fit_a0 <- gam(r1 ~ s(z1) + s(z2) + s(z3) + s(z4) + z5,
  #                   family = binomial,
  #                   data = dat_a0)
  rp1_fit_a1 <- gam(r1 ~ s(x1) + s(x2) + s(x3) + s(x4) + x5,
                    family = binomial,
                    data = dat_a1)
  rp1_fit_a0 <- gam(r1 ~ s(x1) + s(x2) + s(x3) + s(x4) + x5,
                    family = binomial,
                    data = dat_a0)
  pi11_c <- predict(rp1_fit_a1, newdata = dat_trt, type = "response")
  pi10_c <- predict(rp1_fit_a0, newdata = dat_ctl, type = "response")
  
  ## (rp*om ~ x + a) -> Use GAM
  # For R1 = 1, plug in the fitted value of pi2a*mu2a, 
  # regress on z and a
  pimu_value <- function(pi2, mu2){
    pi2*mu2
  }
  
  rpom1_fit_a1 <- gam(pimu_value(pi21_c[dat_R11$a == 1], mu21_c[dat_R11$a == 1]) ~ s(x1) + s(x2) + s(x3) + s(x4) + x5,
                      family = gaussian,
                      data = dat_R11a1, na.action = na.omit)
  rp1nom_fit_a1 <- gam(pimu_value(1 - pi21_c[dat_R11$a == 1], mu20_c[dat_R11$a == 1]) ~ s(x1) + s(x2) + s(x3) + s(x4) + x5,
                       family = gaussian,
                       data = dat_R11a1, na.action = na.omit)
  pimu11_c <- predict(rpom1_fit_a1, newdata = dat_trt)
  pimu1n1_c <- predict(rp1nom_fit_a1, newdata = dat_trt)
  
  rpom_est <- function(pi11, pimu11, pimu1n1, mu10){
    point_value <- mean(pi11*(pimu11 + pimu1n1 - mu10))
    var_value <- var(pi11*(pimu11 + pimu1n1 - mu10))/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  psom_est <- function(e, mu20, mu10){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    r11nr2mu20 <- dat$r1*(1 - dat$r2)
    r11nr2mu20[dat$r1 == 1 & dat$r2 == 0] <- mu20[dat_R11$r2 == 0]
    first_part <- dat$a/e*(r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10)
    second_part <- (1 - dat$a)/(1 - e)*(r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10)
    point_value <- mean(first_part - second_part)
    var_value <- var(first_part - second_part)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  psom_alter_est <- function(e, mu20, mu10){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    r11nr2mu20 <- dat$r1*(1 - dat$r2)
    r11nr2mu20[dat$r1 == 1 & dat$r2 == 0] <- mu20[dat_R11$r2 == 0]
    first_part_up <- dat$a/e*(r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10)
    first_part_down <- dat$a/e
    second_part_up <- (1 - dat$a)/(1 - e)*(r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10)
    second_part_down <- (1 - dat$a)/(1 - e)
    point_value <- sum(first_part_up)/sum(first_part_down) - sum(second_part_up)/sum(second_part_down)
    var_value <- var(first_part_up/mean(first_part_down) - second_part_up/mean(second_part_down))/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  psrp_est <- function(e, pi11, pi10, pi21, pi20, ps_ratio){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    first_part <- dat$a/e*r2y2
    cum_pi20 <- pi10[dat$r1 == 1]*pi20
    prod_part2 <- pi11[dat$r1 == 1]*(1-pi21)*ps_ratio - pi11[dat$r1 == 1]
    frac_part2 <- prod_part2/cum_pi20
    part2 <- dat$r2
    part2[dat$r2 == 1] <- frac_part2[dat_R11$r2 == 1]
    second_part <- (1 - dat$a)/(1 - e)*(part2*r2y2)
    point_value <- mean(first_part + second_part)
    var_value <- var(first_part + second_part)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  psrp_alter_est <- function(e, pi11, pi10, pi21, pi20, ps_ratio){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    first_part_up <- dat$a/e*r2y2
    first_part_down <- dat$a/e
    cum_pi20 <- pi10[dat$r1 == 1]*pi20
    prod_part2 <- pi11[dat$r1 == 1]*(1-pi21)*ps_ratio - pi11[dat$r1 == 1]
    frac_part2 <- prod_part2/cum_pi20
    part2 <- dat$r2
    part2[dat$r2 == 1] <- frac_part2[dat_R11$r2 == 1]
    second_part_up <- (1 - dat$a)/(1 - e)*(part2*r2y2)
    r2_ratio_pi20 <- dat$r2
    r2_ratio_pi20[dat$r2 == 1] <- 1/cum_pi20[dat_R11$r2 == 1]
    second_part_down <- (1 - dat$a)/(1 - e)*r2_ratio_pi20
    point_value <- sum(first_part_up)/sum(first_part_down) + sum(second_part_up)/sum(second_part_down)
    var_value <- var(first_part_up/mean(first_part_down) + second_part_up/mean(second_part_down))/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  psrpom_est <- function(e, mu20, mu10, pimu11, pimu1n1, pi21, pi20, pi11, pi10, ps_ratio){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    r11nr2mu20 <- dat$r1*(1 - dat$r2)
    r11nr2mu20[dat$r1 == 1 & dat$r2 == 0] <- mu20[dat_R11$r2 == 0]
    first_part <- dat$a/e*(r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10)
    second_part <- (1 - dat$a/e)*(pi11*(pimu11 + pimu1n1) + (1 - pi11)*mu10)
    third_part <- mu10
    
    cum_pi20 <- pi10[dat$r1 == 1]*pi20
    prod_part2 <- pi11[dat$r1 == 1]*(1-pi21)*ps_ratio - pi11[dat$r1 == 1]
    frac_part2 <- prod_part2/cum_pi20
    part2 <- dat$r2
    part2[dat$r2 == 1] <- frac_part2[dat_R11$r2 == 1]
    long_mu20 <- dat$r2
    long_mu20[dat$r2 == 1] <- mu20[dat_R11$r2 == 1]
    r2y2_mu20 <- dat$r2*(dat$y2 - long_mu20)
    r2y2_mu20 <- ifelse(is.na(r2y2_mu20), 0, r2y2_mu20)
    
    long1_mu20 <- dat$r1
    long1_mu20[dat$r1 == 1] <- mu20
    r1y1_mu10 <- dat$r1*(long1_mu20 - mu10)
    r1y1_mu10 <- ifelse(is.na(r1y1_mu10), 0, r1y1_mu10)
    
    fourth_part <- (1 - dat$a)/(1 - e)*(r2y2_mu20*part2 - r1y1_mu10*pi11/pi10)
    point_value <- mean(first_part + second_part - third_part + fourth_part)
    var_value <- mean((first_part + second_part - third_part + fourth_part - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  psrpom_alter_est <- function(e, mu20, mu10, pimu11, pimu1n1, pi21, pi20, pi11, pi10, ps_ratio){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    r11nr2mu20 <- dat$r1*(1 - dat$r2)
    r11nr2mu20[dat$r1 == 1 & dat$r2 == 0] <- mu20[dat_R11$r2 == 0]
    first_part_up <- dat$a/e*(r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10)
    first_part_down <- dat$a/e
    second_part1 <- pi11*(pimu11 + pimu1n1) + (1 - pi11)*mu10
    second_part2_up <- dat$a/e*(pi11*(pimu11 + pimu1n1) + (1 - pi11)*mu10)
    second_part2_down <- dat$a/e
    third_part <- mu10
    
    cum_pi20 <- pi10[dat$r1 == 1]*pi20
    prod_part2 <- pi11[dat$r1 == 1]*(1-pi21)*ps_ratio - pi11[dat$r1 == 1]
    frac_part2 <- prod_part2/cum_pi20
    part2 <- dat$r2
    part2[dat$r2 == 1] <- frac_part2[dat_R11$r2 == 1]
    long_mu20 <- dat$r2
    long_mu20[dat$r2 == 1] <- mu20[dat_R11$r2 == 1]
    r2y2_mu20 <- dat$r2*(dat$y2 - long_mu20)
    r2y2_mu20 <- ifelse(is.na(r2y2_mu20), 0, r2y2_mu20)
    
    long1_mu20 <- dat$r1
    long1_mu20[dat$r1 == 1] <- mu20
    r1y1_mu10 <- dat$r1*(long1_mu20 - mu10)
    r1y1_mu10 <- ifelse(is.na(r1y1_mu10), 0, r1y1_mu10)
    
    fourth_part1_up <- (1 - dat$a)/(1 - e)*r2y2_mu20*part2
    
    r2_ratio_pi20 <- dat$r2
    r2_ratio_pi20[dat$r2 == 1] <- 1/cum_pi20[dat_R11$r2 == 1]
    fourth_part1_down <- (1 - dat$a)/(1 - e)*r2_ratio_pi20
    
    fourth_part2_up <- (1 - dat$a)/(1 - e)*r1y1_mu10*pi11/pi10
    fourth_part2_down <- (1 - dat$a)/(1 - e)*dat$r1/pi10
    
    point_value <- sum(first_part_up)/sum(first_part_down) + mean(second_part1) - 
      sum(second_part2_up)/sum(second_part2_down) - mean(third_part) +
      sum(fourth_part1_up)/sum(fourth_part1_down) - sum(fourth_part2_up)/sum(fourth_part2_down)
    var_value <- mean((first_part_up/mean(first_part_down) + second_part1 - 
                         second_part2_up/mean(second_part2_down) - third_part +
                         fourth_part1_up/mean(fourth_part1_down) - fourth_part2_up/mean(fourth_part2_down) - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  ## Calibration
  # (1) Select g(X): include Z, Z^2 and interaction terms
  inter_mat <- matrix(0, nrow = n, ncol = choose(k-1, 2)+k-1)
  count <- 1
  for(i in 1:(k-1)){
    for(j in i:(k-1)){
      inter_mat[,count] <- dat[,k+i]*dat[,k+j]
      count <- count + 1
    }
  }
  g_whole_mat <- cbind(data.matrix(dat[,k + 1:k]), inter_mat)
  tilde_g <- colMeans(g_whole_mat)
  
  # (2) For ps
  g1_mat <- g_whole_mat[which(dat$a == 1),]
  g0_mat <- g_whole_mat[which(dat$a == 0),]
  
  weight_fn <- function(g_mat, tilde_g){
    opt_fn <- function(lambda){
      first_part <- as.vector(exp(g_mat%*%lambda) + 1)/sum(exp(g_mat%*%lambda) + 1)
      colSums(first_part*g_mat) - tilde_g
    }
    set.seed(1234)
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
  
  # For rp
  ## (a) t = 1
  g_whole_mat <- cbind(g_whole_mat, dat$a)
  tilde_ga <- colMeans(g_whole_mat)
  gr1_mat <- g_whole_mat[dat$r1 == 1,]
  wr1_cal <- weight_fn(gr1_mat, tilde_ga)
  
  ## (b) t = 2
  mat_add <- cbind(dat$y1, data.matrix(dat[,k + 1:(k-1)])*dat$y1, dat$y1^2)
  g_whole_mat <- cbind(g_whole_mat, mat_add)
  tilde_gr1 <- colMeans(g_whole_mat[dat$r1 == 1,])
  gr2_mat <- g_whole_mat[dat$r2 == 1,]
  wr2_cal <- weight_fn(gr2_mat, tilde_gr1)
  
  ind1 <- which(dat[which(dat$a == 0), ]$r2 == 1)
  ind2 <- which(dat[which(dat$r1 == 1), ]$a == 0 & dat[which(dat$r1 == 1), ]$r2 == 1)
  ind3 <- which(dat[which(dat$r2 == 1), ]$a == 0)

  w0r2_cal <- (w0_cal[ind1]*wr1_cal[ind2]*wr2_cal[ind3])/sum(w0_cal[ind1]*wr1_cal[ind2]*wr2_cal[ind3])
  
  ind1 <- which(dat[which(dat$a == 0), ]$r1 == 1)
  ind2 <- which(dat[which(dat$r1 == 1), ]$a == 0)
  w0r1_part2_cal <- (w0_cal[ind1]*wr1_cal[ind2])/sum(w0_cal[ind1]*wr1_cal[ind2])
  
  psrpom_weight_est <- function(w1, w0r2, w0r1_part2, 
                                mu20, mu10, pimu11, pimu1n1, pi21,
                                pi11, pi10, ps_ratio){
    r2y2 <- dat$r2*dat$y2
    r2y2 <- ifelse(is.na(r2y2), 0, r2y2)
    r11nr2mu20 <- dat$r1*(1 - dat$r2)
    r11nr2mu20[dat$r1 == 1 & dat$r2 == 0] <- mu20[dat_R11$r2 == 0]
    ind_a1 <- which(dat$a == 1)
    first_part <- r2y2 + r11nr2mu20 + (1 - dat$r1)*mu10
    second_part <- pi11*(pimu11 + pimu1n1) + (1 - pi11)*mu10
    
    sum_part1 <- sum(w1*(first_part[ind_a1] - second_part[ind_a1]))
    sum_part2 <- mean(second_part - mu10)
    
    long_mu20 <- dat$r2
    long_mu20[dat$r2 == 1] <- mu20[dat_R11$r2 == 1]
    r2y2_mu20 <- dat$r2*(dat$y2 - long_mu20)
    r2y2_mu20 <- ifelse(is.na(r2y2_mu20), 0, r2y2_mu20)
    
    third_part <- (ps_ratio[which(dat_R11$r2 == 1 & dat_R11$a == 0)]*pi11[which(dat$r2 == 1 & dat$a == 0)]*(1-pi21)[dat_R11$r2 == 1 & dat_R11$a == 0] - pi11[which(dat$r2 == 1 & dat$a == 0)])*r2y2_mu20[dat$r2 == 1 & dat$a == 0]
    sum_part3 <- sum(w0r2*third_part)
    
    long1_mu20 <- dat$r1
    long1_mu20[dat$r1 == 1] <- mu20
    r1y1_mu10 <- dat$r1*(long1_mu20 - mu10)
    r1y1_mu10 <- ifelse(is.na(r1y1_mu10), 0, r1y1_mu10)
    ind_a0r1 <- which(dat$a == 0 & dat$r1 == 1)
    sum_part4 <- sum(w0r1_part2*(r1y1_mu10*pi11)[ind_a0r1])
    
    point_value <- sum_part1 + sum_part2 + sum_part3 - sum_part4
    
    # if put other weights as 0
    part1_long <- rep(0, n)
    part1_long[dat$a == 1]<- n*w1*(first_part[ind_a1] - second_part[ind_a1])
    part2_long <- second_part - mu10
    part3_long <- rep(0, n)
    part3_long[dat$r2 == 1 & dat$a == 0] <- n*w0r2*third_part
    part4_long <- rep(0, n)
    part4_long[dat$r1 == 1 & dat$a == 0] <- n*w0r1_part2*(r1y1_mu10*pi11)[ind_a0r1]
    
    var_value <- mean((part1_long + part2_long + part3_long - part4_long - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }
  
  # CASE 1: c(ps, om, rp)
  tr_res <- psrpom_est(e_c, mu20_c, mu10_c, pimu11_c,
                       pimu1n1_c, pi21_c, pi20_c,
                       pi11_c, pi10_c, ps_ratio_c)
  tr_psrpom_psrpom_est <- tr_res$point_value
  tr_psrpom_psrpom_var <- tr_res$var_value
  tr_alter_res <- psrpom_alter_est(e_c, mu20_c, mu10_c, pimu11_c,
                                   pimu1n1_c, pi21_c, pi20_c,
                                   pi11_c, pi10_c, ps_ratio_c)
  tr_psrpom_psrpom_alter_est <- tr_alter_res$point_value
  tr_psrpom_psrpom_alter_var <- tr_alter_res$var_value
  psrp_psrpom_res <- psrp_est(e_c, pi11_c, pi10_c, pi21_c, pi20_c,
                              ps_ratio_c)
  psrp_psrpom_est <- psrp_psrpom_res$point_value
  psrp_psrpom_var <- psrp_psrpom_res$var_value
  psrp_psrpom_alter_res <- psrp_alter_est(e_c, pi11_c, pi10_c, pi21_c, pi20_c,
                                          ps_ratio_c)
  psrp_psrpom_alter_est <- psrp_psrpom_alter_res$point_value
  psrp_psrpom_alter_var <- psrp_psrpom_alter_res$var_value
  psom_psrpom_res <- psom_est(e_c, mu20_c, mu10_c)
  psom_psrpom_est <- psom_psrpom_res$point_value
  psom_psrpom_var <- psom_psrpom_res$var_value
  psom_psrpom_alter_res <- psom_alter_est(e_c, mu20_c, mu10_c)
  psom_psrpom_alter_est <- psom_psrpom_alter_res$point_value
  psom_psrpom_alter_var <- psom_psrpom_alter_res$var_value
  rpom_psrpom_res <- rpom_est(pi11_c, pimu11_c, pimu1n1_c, mu10_c)
  rpom_psrpom_est <- rpom_psrpom_res$point_value
  rpom_psrpom_var <- rpom_psrpom_res$var_value
  tr_cal_psrpom_res <- psrpom_weight_est(w1_cal, w0r2_cal, w0r1_part2_cal, 
                                         mu20_c, mu10_c, pimu11_c, pimu1n1_c, pi21_c,
                                         pi11_c, pi10_c, ps_ratio_c)
  tr_cal_psrpom_est <- tr_cal_psrpom_res$point_value
  tr_cal_psrpom_var <- tr_cal_psrpom_res$var_value
  psrpom_all <- c(tr_psrpom_psrpom_est, 
                  tr_psrpom_psrpom_alter_est,
                  tr_cal_psrpom_est,
                  psrp_psrpom_est, 
                  psrp_psrpom_alter_est,
                  psom_psrpom_est, 
                  psom_psrpom_alter_est,
                  rpom_psrpom_est)
  var_all <- c(tr_psrpom_psrpom_var, 
               tr_psrpom_psrpom_alter_var,
               tr_cal_psrpom_var,
               psrp_psrpom_var, 
               psrp_psrpom_alter_var,
               psom_psrpom_var, 
               psom_psrpom_alter_var,
               rpom_psrpom_var)
  
  return(list(psrpom_all = psrpom_all,
              var_all = var_all)
  )
}

nonpara_fn <- function(dat, B, point_est){
  n <- nrow(dat)
  psrpom_all_np <- matrix(0, 8, B)
  psrpom_var <- matrix(0, 8, B)
  t_star <- matrix(0, 8, B)
  for(b in 1:B){
    set.seed(b)
    boot_id <- sample(1:n, size = n, replace = TRUE)
    boot_dat <- dat[boot_id,]
    count_a1r2 <- sum(boot_dat$a == 1 & boot_dat$r2 == 1)
    count_a0r2 <- sum(boot_dat$a == 0 & boot_dat$r2 == 1)
    while(count_a1r2 <= floor(n/10) | count_a0r2 <= floor(n/10)){
      boot_id <- sample(1:n, size = n, replace = TRUE)
      boot_dat <- dat[boot_id,]
      count_a1r2 <- sum(boot_dat$a == 1 & boot_dat$r2 == 1)
      count_a0r2 <- sum(boot_dat$a == 0 & boot_dat$r2 == 1)
    }
    boot_res <- model_est(boot_dat)
    psrpom_all_np[,b] <- boot_res$psrpom_all
    psrpom_var[,b] <- boot_res$var_all
    t_star[,b] <- (psrpom_all_np[,b] - point_est)/sqrt(psrpom_var[,b])
  }
  c_star <- apply(abs(t_star), 1, function(x) quantile(x, 0.95))
  var_boot <- apply(psrpom_all_np, 1, var)
  percentile_boot <- apply(psrpom_all_np, 1, function(x) quantile(x, c(0.025, 0.975)))
  
  return(list(c_star = c_star,
              var_boot = var_boot,
              percentile_boot = percentile_boot))
}

main <- function(seed){
  set.seed(seed)
  dat_list <- sim_dat(n = n, k = k, alpha = alpha, gamma1 = gamma1, 
                      gamma2 = gamma2)
  dat <- dat_list$dat
  true_value <- dat_list$true_value
  true_value_psom <- dat_list$true_value_psom
  psrpom_res <- model_est(dat)
  psrpom_all <- psrpom_res$psrpom_all
  psrpom_var <- psrpom_res$var_all
  boot_res <- nonpara_fn(dat, B, psrpom_all)
  c_star <- boot_res$c_star
  var_boot <- boot_res$var_boot
  percentile_boot <- boot_res$percentile_boot
  
  return(list(point_est = psrpom_all,
              var_est = psrpom_var,
              # true value
              true_value = true_value_psom,
              # c-star
              c_star = c_star,
              var_boot = var_boot,
              percentile_boot = percentile_boot
  ))
}
