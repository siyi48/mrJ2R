#' Obtain the proposed estimator in the cross-sectional study
#'
#' This function derives the estimators for the average treatment effect (ATE)
#' under jump-to-reference (J2R) in the cross-sectional study. The resulting ATE
#' estimators consists of five simple estimators (rp-om, ps-om, ps-om-N, ps-rp,
#' and ps-rp-N) and three triply robust estimators (tr, tr-N, and tr-C)
#' motivated by the efficient influence function (EIF). Also, EIF-based
#' variance estimators are provided for the three EIF-based estimators. The
#' nuisance functions are estimated by parametric models.
#'
#' @param formula.ps a list of formulas for fitting the propensity score model
#' @param formula.rp a list of formulas for fitting the response probability
#' model
#' @param formula.om a list of formulas for fitting the outcome mean model
#' @param mat.cal the matrix for calibration, should have the same subject order
#' as the original data.
#' @param data the data that include the binary treatment assignment, where the
#' active treatment group should be encoded as `1`, and the control group should be
#' encoded as `0`, the binary response status, where `1` indicates the subject
#' is observed and `0` indicates the subject is missing, the outcome, and the
#' fully observed baseline covariates.
#' @param type the type of estimate the user want to obtain. Available types
#' include:
#' \itemize{
#' \item `tr`: the triply robust estimator
#' \item `tr.norm`: the triply robust estimator with normalization
#' \item `tr.cal`: the triply robust estimator with calibration
#' \item `psom`: the ps-om estimator
#' \item `psom.norm`: the ps-om estimator with normalization
#' \item `psrp`: the ps-rp estimator
#' \item `psrp.norm`: the ps-rp estimator with normalization
#' \item `rpom`: the rp-om estimator
#' }
#'
#' @return The point estimate
#' @import nleqslv
#' @export
#'
#' @examples
#' data <- dat1
#' formula.om <- y ~ z1 + z2 + z3 + z4 + z5 + a +
#'   a:(z1 + z2 + z3 + z4 + z5)
#' formula.ps <- a ~ z1 + z2 + z3 + z4 + z5
#' formula.rp <- r ~ z1 + z2 + z3 + z4 + z5 + a +
#'   a:(z1 + z2 + z3 + z4 + z5)
#' mat.cal <- data.matrix(data[,paste0("z", 1:5)]) # calibration matrix
#' res.1time <- jtr.1time(formula.ps = formula.ps,
#'                        formula.rp = formula.rp,
#'                        formula.om = formula.om,
#'                        data = data, mat.cal = mat.cal,
#'                        type = c("tr", "tr.cal", "rpom"))
#' res.1time
jtr.1time <- function(formula.ps, formula.rp, formula.om, data,
                      mat.cal = NULL, type){
  # get the data for prediction
  trt.name <- as.character(formula.ps[[2]])
  rp.name <- as.character(formula.rp[[2]])
  y.name <- as.character(formula.om[[2]])
  dat.ctl <- data
  dat.ctl[,trt.name] <- 0
  dat.trt <- data
  dat.trt[,trt.name] <- 1
  a <- data[,trt.name]
  r <- data[,rp.name]
  y <- data[,y.name]
  n <- nrow(data)

  # fit om
  om.fit <- lm(formula.om, data = data, na.action = na.omit)

  mu1 <- predict(om.fit, newdata = dat.trt)
  mu0 <- predict(om.fit, newdata = dat.ctl)
  # fit ps
  ps.fit <- glm(formula.ps, family = binomial, data = data)
  e <- predict(ps.fit, newdata = data, type = "response")
  # fit rp
  rp.fit <- glm(formula.rp, family = binomial, data = data)
  pi1 <- predict(rp.fit, newdata = dat.trt, type = "response")
  pi0 <- predict(rp.fit, newdata = dat.ctl, type = "response")

  rpom_est <- function(pi1, mu1, mu0){
    return(mean(pi1*(mu1 - mu0)))
  }

  psom_est <- function(e, mu1, mu0){
    ry <- r*y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part <- a/e*(ry + (1 - r)*mu0)
    second_part <- (1 - a)/(1 - e)*(ry + (1 - r)*mu0)
    return(mean(first_part - second_part))
  }

  psom_alter_est <- function(e, mu1, mu0){
    ry <- r*y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part_up <- a/e*(ry + (1 - r)*mu0)
    first_part_down <- a/e
    second_part_up <- (1 - a)/(1 - e)*(ry + (1 - r)*mu0)
    second_part_down <- (1 - a)/(1 - e)
    return(sum(first_part_up)/sum(first_part_down) - sum(second_part_up)/sum(second_part_down))
  }

  psrp_est <- function(e, pi1, pi0){
    ry <- r*y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part <- a/e*ry
    second_part <- (1 - a)/(1 - e)*(pi1/pi0*ry)
    return(mean(first_part - second_part))
  }

  psrp_alter_est <- function(e, pi1, pi0){
    ry <- r*y
    ry <- ifelse(is.na(ry), 0, ry)
    first_part_up <- a/e*ry
    first_part_down <- a/e
    second_part_up <- (1 - a)/(1 - e)*(pi1/pi0*ry)
    second_part_down <- (1 - a)*r/((1 - e)*pi0)
    return(sum(first_part_up)/sum(first_part_down) - sum(second_part_up)/sum(second_part_down))
  }

  psrpom_est <- function(e, pi1, pi0, mu1, mu0){
    ry_mu0 <- r*(y - mu0)
    ry_mu0 <- ifelse(is.na(ry_mu0), 0, ry_mu0)
    first_part <- (a/e - (1 - a)/(1 - e)*pi1/pi0)*ry_mu0
    second_part <- (a - e)/e*pi1*(mu1 - mu0)
    point_value <- mean(first_part - second_part)
    var_value <- mean((first_part - second_part - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }

  psrpom_alter_est <- function(e, pi1, pi0, mu1, mu0){
    ry_mu0 <- r*(y - mu0)
    ry_mu0 <- ifelse(is.na(ry_mu0), 0, ry_mu0)
    first_part_up <- a/e*(ry_mu0 - pi1*(mu1 - mu0))
    first_part_down <- a/e
    second_part_up <- (1 - a)/(1 - e)*pi1/pi0*ry_mu0
    second_part_down <- (1 - a)*r/((1 - e)*pi0)
    third_part <- pi1*(mu1 - mu0)
    point_value <- sum(first_part_up)/sum(first_part_down) -
      sum(second_part_up)/sum(second_part_down) + mean(third_part)
    var_value <- mean((first_part_up/mean(first_part_down) -
                         second_part_up/mean(second_part_down) + third_part - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }

  psrpom_weight_est <- function(w1, w0r, pi1, mu1, mu0){
    ry_mu0 <- r*(y - mu0)
    ry_mu0 <- ifelse(is.na(ry_mu0), 0, ry_mu0)
    diff <- ry_mu0 - pi1*(mu1 - mu0)
    diff1 <- diff[which(a == 1)]
    piry_mu00r <- (pi1*ry_mu0)[which(a == 0 & r == 1)]
    first_part <- sum(w1*diff1)
    second_part <- sum(w0r*piry_mu00r)
    third_part <- mean(pi1*(mu1 - mu0))
    point_value <- first_part - second_part + third_part

    part1_long <- rep(0, n)
    part1_long[a == 1]<- n*w1*diff1
    part2_long <- rep(0, n)
    part2_long[a == 0 & r == 1]<- n*w0r*piry_mu00r
    part3_long <- pi1*(mu1 - mu0)

    var_value <- mean((part1_long - part2_long + part3_long - point_value)^2)/n
    return(list(point_value = point_value,
                var_value = var_value))
  }

  est <- NULL
  ve <- NULL
  if("tr" %in% type){
    res <- psrpom_est(e, pi1, pi0, mu1, mu0)
    est <- c(est, res$point_value)
    ve <- c(ve, res$var_value)
    names(ve)[length(ve)] <- "tr"
  }

  if("tr.norm" %in% type){
    res <- psrpom_alter_est(e, pi1, pi0, mu1, mu0)
    est <- c(est, res$point_value)
    ve <- c(ve, res$var_value)
    names(ve)[length(ve)] <- "tr.norm"
  }

  if("tr.cal" %in% type){
    if(is.null(mat.cal)){
      stop("For calibration-based estimators, need to specify `mat.cal`")
    }
    g_whole_mat <- data.matrix(mat.cal)
    tilde_g <- colMeans(g_whole_mat)

    g1_mat <- g_whole_mat[a == 1,]
    g0_mat <- g_whole_mat[a == 0,]

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

    gr_mat <- g_whole_mat[r == 1,]
    wr_cal <- weight_fn(gr_mat, tilde_g)
    temp.mat <- data.frame(a = a, r = r)
    ind_temp1 <- which(temp.mat[which(temp.mat$a == 0), ]$r == 1)
    ind_temp2 <- which(temp.mat[which(temp.mat$r == 1), ]$a == 0)
    w0r_cal <- (w0_cal[ind_temp1]*wr_cal[ind_temp2])/sum(w0_cal[ind_temp1]*wr_cal[ind_temp2])

    res <- psrpom_weight_est(w1_cal, w0r_cal, pi1, mu1, mu0)
    est <- c(est, res$point_value)
    ve <- c(ve, res$var_value)
    names(ve)[length(ve)] <- "tr.cal"
  }

  if("rpom" %in% type){
    est <- c(est, rpom_est(pi1, mu1, mu0))
  }

  if("psom" %in% type){
    est <- c(est, psom_est(e, mu1, mu0))
  }

  if("psom.norm" %in% type){
    est <- c(est, psom_alter_est(e, mu1, mu0))
  }

  if("psrp" %in% type){
    est <- c(est, psrp_est(e, pi1, pi0))
  }

  if("psom.norm" %in% type){
    est <- c(est, psrp_alter_est(e, pi1, pi0))
  }
  names(est) <- type

  return(list(est = est,
              ve = ve))
}


#' Obtain the proposed estimator in the longitudinal study
#'
#' This function derives the estimators for the average treatment effect (ATE)
#' under jump-to-reference (J2R) in the longitudinal study. The resulting ATE
#' estimators consists of five simple estimators (rp-pm, ps-om, ps-om-N, ps-rp,
#' and ps-rp-N) and two multiply robust estimators (mr and mr-N) motivated
#' by the efficient influence function (EIF). Also, EIF-based variance
#' estimators are provided for the two EIF-based estimators. The
#' nuisance functions are estimated by generalized additive models (GAMs).
#' Users can select the type of smoothness used in GAMs.
#'
#' @param formula.ps a list of formulas for fitting the propensity score model.
#' @param formula.rp a list of formulas for fitting the response probability
#' model.
#' @param formula.om a list of formulas for fitting the outcome mean model.
#' @param formula.pm a list of formulas for fitting the pattern mean model.
#' Users can input any variables in the data as the response variables, since
#' the response variables in the pattern mean formulas are recalculated by the
#' estimated outcome mean and response probability.
#' @param list.cal a list of matrices for calibration. Each matrix in `list.cal`
#' should have the same number of rows with the same subject orders. The first
#' matrix in the list should include the variables that may affect the treatment
#' assignment, the latter matrices in the list should include the variables that
#' may affect the response status in a time order.
#' @param data the data that include the binary treatment assignment, where the
#' active treatment group should be encoded as `1`, and the control group should be
#' encoded as `0`, the binary response status, where `1` indicates the subject
#' is observed and `0` indicates the subject is missing, the outcome, and the
#' fully observed baseline covariates.
#' @param type the type of estimate the user want to obtain. Available types
#' include:
#' \itemize{
#' \item `mr`: the multiply robust estimator
#' \item `mr.norm`: the multiply robust estimator with normalization
#' \item `mr.cal`: the multiply robust estimator with calibration
#' \item `psom`: the ps-om estimator
#' \item `psom.norm`: the ps-om estimator with normalization
#' \item `psrp`: the ps-rp estimator
#' \item `psrp.norm`: the ps-rp estimator with normalization
#' \item `rppm`: the rp-pm estimator
#' }
#'
#' @return The point estimate
#' @import nleqslv gam
#' @export
#'
#' @examples
#' data <- dat2
#' formula.om <- list()
#' formula.om[[2]] <- y2 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
#' x5 + gam::s(y1)
#' formula.om[[1]] <- y1 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5
#' formula.ps <- list()
#' formula.ps[[2]] <- a ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
#' x5  + gam::s(y1)
#' formula.ps[[1]] <- a ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5
#' formula.rp <- list()
#' formula.rp[[2]] <- r2 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
#' x5 + gam::s(y1)
#' formula.rp[[1]] <- r1 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5
#' formula.pm <- list()
#' formula.pm[[2]] <- y2 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) +
#' x5 + gam::s(y1)
#' formula.pm[[1]] <- y1 ~ gam::s(x1) + gam::s(x2) + gam::s(x3) + gam::s(x4) + x5
#' k <- 5
#' inter_mat <- matrix(0, nrow = nrow(data), ncol = choose(k-1, 2)+k-1)
#' count <- 1
#' for(i in 1:(k-1)){
#'   for(j in i:(k-1)){
#'     inter_mat[,count] <- data[,k+i]*data[,k+j]
#'     count <- count + 1
#'   }
#' }
#' list.cal <- list()
#' list.cal[[1]] <- cbind(data.matrix(data[,k + 1:k]), inter_mat) # A
#' list.cal[[2]] <- cbind(cbind(list.cal[[1]], data$a)) # R1
#' list.cal[[3]] <- cbind(list.cal[[2]], data$y1,
#'                        data.matrix(data[,k + 1:(k-1)])*data$y1, data$y1^2) # R2
#' res.longi <- jtr.longi(formula.ps = formula.ps,
#'                        formula.rp = formula.rp,
#'                        formula.om = formula.om,
#'                        formula.pm = formula.om,
#'                        list.cal = list.cal,
#'                        data = data,
#'                        type = c("mr", "mr.cal", "rppm"))
#' res.longi
jtr.longi <- function(formula.ps, formula.rp, formula.om, formula.pm,
                      list.cal = NULL,
                      data, type){
  t.time <- length(formula.om)
  trt.name <- as.character(formula.ps[[1]][2])
  rp.name <- as.character(sapply(formula.rp, function(x) x[[2]]))
  y.name <- as.character(sapply(formula.om, function(x) x[[2]]))
  pm.name <- as.character(sapply(formula.pm, function(x) x[[2]]))
  a <- data[,trt.name]
  r.mat <- cbind(1, sapply(rp.name, function(x) data[,x]))
  colnames(r.mat) <- paste0("r", 0:t.time)
  y.mat <- sapply(y.name, function(x) data[,x])
  colnames(y.mat) <- paste0("y", 1:t.time)
  n <- nrow(data)

  # get the data for prediction
  # observed data at the second last visit
  dat.ctl <- data
  dat.ctl[,trt.name] <- 0
  dat.trt <- data
  dat.trt[,trt.name] <- 1
  dat.a1 <- data[a == 1,]
  dat.a0 <- data[a == 0,]

  # fit om
  om.fit.a1 <- gam::gam(formula.om[[t.time]],
                        data = dat.a1, na.action = na.omit)
  dat.obs <- data[r.mat[,t.time] == 1,]
  dat.obsa1 <- dat.obs[dat.obs[,trt.name] == 1,]
  mu.last.a1 <- predict(om.fit.a1, newdata = dat.obs)
  dat.obsa0 <- dat.obs[dat.obs[,trt.name] == 0,]
  om.fit.a0 <- list()
  mu.a0 <- list()
  mu.a0.total <- list()
  for(k in t.time:1){
    dat.obs <- data[r.mat[,k+1] == 1,]
    dat.obsa0 <- dat.obs[dat.obs[,trt.name] == 0,]
    dat.obs.pred <- data[r.mat[,k] == 1,]
    dat.obsa0.pred <- dat.obs.pred[dat.obs.pred[,trt.name] == 0,]
    if(k < t.time){
      dat.obsa0[,y.name[k]] <- mu.a0[[t.time - k]]
    }
    om.fit.a0[[t.time - k + 1]] <- gam::gam(formula.om[[k]],
                                            data = dat.obsa0,
                                            na.action = na.omit)
    mu.a0[[t.time - k + 1]] <- predict(om.fit.a0[[t.time - k + 1]],
                                       newdata = dat.obsa0.pred)
    mu.a0.total[[t.time - k + 1]] <- predict(om.fit.a0[[t.time - k + 1]],
                                             newdata = dat.obs.pred)
  }

  # fit ps
  dat.obs <- data
  ps.fit <- list()
  e.ps <- list()
  ps.ratio <- list()
  ps.ratio[[1]] <- rep(1, n)
  for(k in 1:t.time){
    dat.obs <- data[r.mat[,k] == 1,]
    ps.fit[[k]] <- gam::gam(formula.ps[[k]],
                            family = binomial,
                            data = dat.obs,
                            na.action = na.omit)
    e.ps[[k]] <- predict(ps.fit[[k]], newdata = dat.obs,
                         type = "response")
    if(k > 1){
      ps.ratio[[k]] <- (e.ps[[k]]/(e.ps[[1]][r.mat[,k] == 1]))/((1 - e.ps[[k]])/(1 - e.ps[[1]][r.mat[,k] == 1]))
    }
  }

  # fit rp
  dat.obs <- data[r.mat[,t.time] == 1,]
  dat.obsa1 <- dat.obs[dat.obs[,trt.name] == 1,]
  dat.obsa0 <- dat.obs[dat.obs[,trt.name] == 0,]
  rp.fit.a1 <- list()
  rp.fit.a0 <- list()
  pi.a1 <- list() # for individuals with R_{s-1} = 1 and A = 1
  pi.a1.total <- list() # for all individuals with R_{s-1} = 1
  pi.a0 <- list() # for individuals with R_{s-1} = 1 and A = 0
  pi.a0.total <- list() # for all individuals with R_{s-1} = 1 and A = 0
  for(k in t.time:1){
    dat.obs <- data[r.mat[,k] == 1,]
    dat.obsa1 <- dat.obs[dat.obs[,trt.name] == 1,]
    dat.obsa0 <- dat.obs[dat.obs[,trt.name] == 0,]
    # For those with A = 1
    rp.fit.a1[[t.time - k + 1]] <- gam::gam(formula.rp[[k]],
                                            data = dat.obsa1,
                                            family = binomial,
                                            na.action = na.omit)
    pi.a1[[t.time - k + 1]] <- predict(rp.fit.a1[[t.time - k + 1]],
                                       newdata = dat.obsa1,
                                       type = "response")
    pi.a1.total[[t.time - k + 1]] <- predict(rp.fit.a1[[t.time - k + 1]],
                                             newdata = dat.obs,
                                             type = "response")
    # For those with A = 0
    rp.fit.a0[[t.time - k + 1]] <- gam::gam(formula.rp[[k]],
                                            data = dat.obsa0,
                                            family = binomial,
                                            na.action = na.omit)
    pi.a0[[t.time - k + 1]] <- predict(rp.fit.a0[[t.time - k + 1]],
                                       newdata = dat.obsa0,
                                       type = "response")
    pi.a0.total[[t.time - k + 1]] <- predict(rp.fit.a0[[t.time - k + 1]],
                                             newdata = dat.obs,
                                             type = "response")
  }

  # fit pm
  pimu_value <- function(pi, mu){
    pi*mu
  }

  pm.a1 <- list()
  for(k in t.time:1){
    dat.obs <- data[r.mat[,k+1] == 1,]
    dat.obsa1 <- dat.obs[dat.obs[,trt.name] == 1,]
    dat.obs.pred <- data[r.mat[,k] == 1,]
    pm.a1[[t.time - k + 1]] <- matrix(0, nrow = nrow(dat.obs.pred),
                                      ncol = t.time - k + 1)
    for(j in 1:(t.time - k + 1)){
      if(k < t.time & j < (t.time - k + 1)){
        dat.obsa1[,pm.name[k]] <- pimu_value(pi.a1[[t.time - k]],
                                            pm.a1[[t.time - k]][dat.obs[,trt.name] == 1,j])
      }
      else if(k < t.time & j == (t.time - k + 1)){
        mu0.a1 <- mu.a0.total[[j-1]][dat.obs[,trt.name] == 1]
        dat.obsa1[,pm.name[k]] <- pimu_value((1 - pi.a1[[t.time - k]]), mu0.a1)
      }
      pm.fit <- gam::gam(formula.pm[[k]],
                         data = dat.obsa1,
                         na.action = na.omit)
      pm.a1[[t.time - k + 1]][,j] <- predict(pm.fit, newdata = dat.obs.pred)
    }
  }

  # rppm_est(pi.a1.total, pm.a1, mu.a0.total)
  rppm_est <- function(pi.a1.total, pm.a1, mu.a0.total){
    point_value <- mean(pi.a1.total[[t.time]]*(rowSums(pm.a1[[t.time]]) -
                                                 mu.a0.total[[t.time]]))
    return(point_value)
  }

  # psom_est(e.ps, mu.a0.total)
  psom_est <- function(e.ps, mu.a0.total){

    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp)
    temp.value <- rep(0, n)
    for(k in t.time:1){
      y.temp <- rep(0, n)
      y.temp[r.mat[,t.time - k + 1] == 1] <- mu.a0.total[[k]]

      temp.value <- r.mat[,t.time - k + 1]*(1 - r.mat[,t.time - k + 2])*y.temp + temp.value
    }
    y.imp <- y.imp + temp.value

    first_part <- a/e.ps[[1]]*y.imp
    second_part <- (1 - a)/(1 - e.ps[[1]])*y.imp
    point_value <- mean(first_part - second_part)
    return(point_value)
  }

  # psom_alter_est(e.ps, mu.a0.total)
  psom_alter_est <- function(e.ps, mu.a0.total){
    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp)
    temp.value <- rep(0, n)
    for(k in t.time:1){
      y.temp <- rep(0, n)
      y.temp[r.mat[,t.time - k + 1] == 1] <- mu.a0.total[[k]]

      temp.value <- r.mat[,t.time - k + 1]*(1 - r.mat[,t.time - k + 2])*y.temp + temp.value
    }
    y.imp <- y.imp + temp.value

    first_part_up <- a/e.ps[[1]]*y.imp
    first_part_down <- a/e.ps[[1]]
    second_part_up <- (1 - a)/(1 - e.ps[[1]])*y.imp
    second_part_down <- (1 - a)/(1 - e.ps[[1]])
    point_value <- sum(first_part_up)/sum(first_part_down) - sum(second_part_up)/sum(second_part_down)
    return(point_value)
  }

  # psrp_est(e.ps, pi.a0.total, pi.a1.total, ps.ratio)
  psrp_est <- function(e.ps, pi.a0.total, pi.a1.total, ps.ratio){
    n.complete <- sum(r.mat[,t.time+1])
    index.complete <- which(r.mat[,t.time+1] == 1)
    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp) # RtYt

    pi.a0.complete <- matrix(0, nrow = n.complete, ncol = t.time)
    pi.a1.complete <- matrix(0, nrow = n.complete, ncol = t.time)
    ps.ratio.complete <- matrix(0, nrow = n.complete, ncol = t.time)
    for(k in t.time:1){
      index <- which(r.mat[r.mat[,k] == 1, t.time+1] == 1)
      pi.a0.complete[,k] <- pi.a0.total[[t.time - k + 1]][index]
      pi.a1.complete[,k] <- pi.a1.total[[t.time - k + 1]][index]
      ps.ratio.complete[,k] <- ps.ratio[[k]][index]
    }
    pi.a0.cum <- cbind(1, t(apply(pi.a0.complete, 1, cumprod)))
    pi.a0.cum.long <- rep(1, n)
    pi.a0.cum.long[index.complete] <- pi.a0.cum[,t.time+1]

    sum.part2 <- rep(0, n.complete)
    for(k in t.time:1){
      sum.part2 <- pi.a0.cum[,k]*(1 - pi.a1.complete[,k])*ps.ratio.complete[,k] + sum.part2
    }
    sum.part2.long <- rep(0, n)
    sum.part2.long[index.complete] <- sum.part2

    part2 <- (sum.part2.long - 1)*(1 - a)*y.imp/((1 - e.ps[[1]])*pi.a0.cum.long)
    part1 <- a*y.imp/e.ps[[1]]

    point_value <- mean(part1 + part2)
    return(point_value)
  }

  # psrp_alter_est(e.ps, pi.a0.total, pi.a1.total, ps.ratio)
  psrp_alter_est <- function(e.ps, pi.a0.total, pi.a1.total, ps.ratio){
    n.complete <- sum(r.mat[,t.time+1])
    index.complete <- which(r.mat[,t.time+1] == 1)
    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp) # RtYt

    pi.a0.complete <- matrix(0, nrow = n.complete, ncol = t.time)
    pi.a1.complete <- matrix(0, nrow = n.complete, ncol = t.time)
    ps.ratio.complete <- matrix(0, nrow = n.complete, ncol = t.time)
    for(k in t.time:1){
      index <- which(r.mat[r.mat[,k] == 1, t.time+1] == 1)
      pi.a0.complete[,k] <- pi.a0.total[[t.time - k + 1]][index]
      pi.a1.complete[,k] <- pi.a1.total[[t.time - k + 1]][index]
      ps.ratio.complete[,k] <- ps.ratio[[k]][index]
    }
    pi.a0.cum <- cbind(1, t(apply(pi.a0.complete, 1, cumprod)))
    pi.a0.cum.long <- rep(1, n)
    pi.a0.cum.long[index.complete] <- pi.a0.cum[,t.time+1]

    sum.part2 <- rep(0, n.complete)
    for(k in t.time:1){
      sum.part2 <- pi.a0.cum[,k]*(1 - pi.a1.complete[,k])*ps.ratio.complete[,k] + sum.part2
    }
    sum.part2.long <- rep(0, n)
    sum.part2.long[index.complete] <- sum.part2

    part2 <- (sum.part2.long - 1)*(1 - a)*y.imp/((1 - e.ps[[1]])*pi.a0.cum.long)
    part1 <- a*y.imp/e.ps[[1]]

    first_part_up <- part1
    first_part_down <- a/e.ps[[1]]

    second_part_up <- part2
    second_part_down <- (1 - a)*r.mat[,ncol(r.mat)]/((1 - e.ps[[1]])*pi.a0.cum.long)

    return(sum(first_part_up)/sum(first_part_down) + sum(second_part_up)/sum(second_part_down))
  }

  # mr_est(e.ps, pm.a1, mu.a0.total,
  #        pi.a0.total, pi.a1.total, ps.ratio)
  mr_est <- function(e.ps, pm.a1, mu.a0.total,
                     pi.a0.total, pi.a1.total, ps.ratio){
    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp)
    temp.value <- rep(0, n)
    for(k in t.time:1){
      y.temp <- rep(0, n)
      y.temp[r.mat[,t.time - k + 1] == 1] <- mu.a0.total[[k]]

      temp.value <- r.mat[,t.time - k + 1]*(1 - r.mat[,t.time - k + 2])*y.temp + temp.value
    }
    y.imp <- y.imp + temp.value

    first_part <- a/e.ps[[1]]*y.imp

    second_part <- (1 - a/e.ps[[1]])*(pi.a1.total[[t.time]]*(rowSums(pm.a1[[t.time]]) -
                                                             mu.a0.total[[t.time]]) + mu.a0.total[[t.time]])
    third_part <- mu.a0.total[[t.time]]

    n.complete <- sum(r.mat[,t.time+1])
    index.complete <- which(r.mat[,t.time+1] == 1)

    pi.a0.obs <- matrix(0, nrow = n, ncol = t.time)
    pi.a1.obs <- matrix(0, nrow = n, ncol = t.time)
    ps.ratio.obs <- matrix(0, nrow = n, ncol = t.time)
    mu.a0.obs <- matrix(0, nrow = n, ncol = t.time+1)
    for(k in t.time:1){
      index <- which(r.mat[,k] == 1)
      pi.a0.obs[index,k] <- pi.a0.total[[t.time - k + 1]]
      pi.a1.obs[index,k] <- pi.a1.total[[t.time - k + 1]]
      ps.ratio.obs[index,k] <- ps.ratio[[k]]
      mu.a0.obs[index,k] <- mu.a0.total[[t.time - k + 1]]
    }
    pi.a0.cum <- cbind(1, t(apply(pi.a0.obs, 1, cumprod)))
    pi.a0.cum <- ifelse(pi.a0.cum == 0, 2, pi.a0.cum)
    mu.a0.obs[index.complete, t.time + 1] <- data[index.complete,y.name[t.time]]
    pi.a1.obs <- cbind(1, pi.a1.obs)

    # recode: cumulative probability (using all observed),
    # sum.part2 / sum.part.seq
    sum.part2 <- rep(0, n)
    sum.part.seq <- list()
    sum.part.seq.long <- matrix(0, nrow = n, ncol = t.time)
    mu.diff <- matrix(0, nrow = n, ncol = t.time)
    for(k in 1:t.time){
      sum.part2 <- pi.a0.cum[,k]*(1 - pi.a1.obs[,k+1])*ps.ratio.obs[,k] + sum.part2
      sum.part.seq[[k]] <- sum.part2
      mu.diff[,k] <- mu.a0.obs[,t.time - k + 2] - mu.a0.obs[,t.time - k + 1]
    }

    sum.part4 <- rep(0, n)
    for(k in 1:t.time){
      index <- which(r.mat[,k] == 1)
      sum.part4 <- (sum.part.seq[[k]] - 1)*mu.diff[,t.time - k + 1]*r.mat[,k+1]/pi.a0.cum[,k+1] + sum.part4
    }

    fourth_part <- (1 - a)*sum.part4/(1 - e.ps[[1]])

    point_value <- mean(first_part + second_part - third_part + fourth_part)
    var_value <- mean((first_part + second_part - third_part + fourth_part - point_value)^2)/n

    return(list(point_value = point_value,
                var_value = var_value))
  }

  # mr_alter_est(e.ps, pm.a1, mu.a0.total,
  #        pi.a0.total, pi.a1.total, ps.ratio)
  mr_alter_est <- function(e.ps, pm.a1, mu.a0.total,
                           pi.a0.total, pi.a1.total, ps.ratio){
    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp)
    temp.value <- rep(0, n)
    for(k in t.time:1){
      y.temp <- rep(0, n)
      y.temp[r.mat[,t.time - k + 1] == 1] <- mu.a0.total[[k]]

      temp.value <- r.mat[,t.time - k + 1]*(1 - r.mat[,t.time - k + 2])*y.temp + temp.value
    }
    y.imp <- y.imp + temp.value

    first_part_up <- a/e.ps[[1]]*y.imp
    first_part_down <- a/e.ps[[1]]

    second_part1 <- pi.a1.total[[t.time]]*(rowSums(pm.a1[[t.time]]) -
                                             mu.a0.total[[t.time]]) + mu.a0.total[[t.time]]
    second_part2_up <- a/e.ps[[1]]*second_part1
    second_part2_down <- a/e.ps[[1]]

    third_part <- mu.a0.total[[t.time]]

    n.complete <- sum(r.mat[,t.time+1])
    index.complete <- which(r.mat[,t.time+1] == 1)

    pi.a0.obs <- matrix(0, nrow = n, ncol = t.time)
    pi.a1.obs <- matrix(0, nrow = n, ncol = t.time)
    ps.ratio.obs <- matrix(0, nrow = n, ncol = t.time)
    mu.a0.obs <- matrix(0, nrow = n, ncol = t.time+1)
    for(k in t.time:1){
      index <- which(r.mat[,k] == 1)
      pi.a0.obs[index,k] <- pi.a0.total[[t.time - k + 1]]
      pi.a1.obs[index,k] <- pi.a1.total[[t.time - k + 1]]
      ps.ratio.obs[index,k] <- ps.ratio[[k]]
      mu.a0.obs[index,k] <- mu.a0.total[[t.time - k + 1]]
    }
    pi.a0.cum <- cbind(1, t(apply(pi.a0.obs, 1, cumprod)))
    pi.a0.cum <- ifelse(pi.a0.cum == 0, 2, pi.a0.cum)
    mu.a0.obs[index.complete, t.time + 1] <- data[index.complete,y.name[t.time]]
    pi.a1.obs <- cbind(1, pi.a1.obs)

    # recode: cumulative probability (using all observed),
    # sum.part2 / sum.part.seq
    sum.part2 <- rep(0, n)
    sum.part.seq <- list()
    sum.part.seq.long <- matrix(0, nrow = n, ncol = t.time)
    mu.diff <- matrix(0, nrow = n, ncol = t.time)
    for(k in 1:t.time){
      sum.part2 <- pi.a0.cum[,k]*(1 - pi.a1.obs[,k+1])*ps.ratio.obs[,k] + sum.part2
      sum.part.seq[[k]] <- sum.part2
      mu.diff[,k] <- mu.a0.obs[,t.time - k + 2] - mu.a0.obs[,t.time - k + 1]
    }

    sum.part4 <- rep(0, n)
    mean.part4 <- rep(0, n)
    fourth_part_up <- matrix(0, nrow = n, ncol = t.time)
    fourth_part_down <- matrix(0, nrow = n, ncol = t.time)
    for(k in 1:t.time){
      fourth_part_up <- (1 - a)*(sum.part.seq[[k]] - 1)*mu.diff[,t.time - k + 1]*r.mat[,k+1]/((1 - e.ps[[1]])*pi.a0.cum[,k+1])
      fourth_part_down <- (1 - a)*r.mat[,k+1]/((1 - e.ps[[1]])*pi.a0.cum[,k+1])
      sum.part4 <- fourth_part_up/sum(fourth_part_down) + sum.part4
      mean.part4 <- fourth_part_up/mean(fourth_part_down) + mean.part4
    }

    fourth_part <- sum.part4

    point_value <- sum(first_part_up)/sum(first_part_down) + mean(second_part1) -
      sum(second_part2_up)/sum(second_part2_down) - mean(third_part) + sum(fourth_part)
    var_value <- mean((first_part_up/mean(first_part_down) + second_part1 -
                         second_part2_up/mean(second_part2_down) - third_part +
                         mean.part4 - point_value)^2)/n

    return(list(point_value = point_value,
                var_value = var_value))
  }

  # mr_weight_est(w.cal, pm.a1, mu.a0.total,pi.a0.total, pi.a1.total, ps.ratio)
  mr_weight_est <- function(w.cal, pm.a1, mu.a0.total,
                            pi.a0.total, pi.a1.total, ps.ratio){
    y.imp <- r.mat[,ncol(r.mat)]*data[,y.name[t.time]]
    y.imp <- ifelse(is.na(y.imp), 0, y.imp)
    temp.value <- rep(0, n)
    for(k in t.time:1){
      y.temp <- rep(0, n)
      y.temp[r.mat[,t.time - k + 1] == 1] <- mu.a0.total[[k]]
      temp.value <- r.mat[,t.time - k + 1]*(1 - r.mat[,t.time - k + 2])*y.temp + temp.value
    }
    y.imp <- y.imp + temp.value

    second_part1 <- pi.a1.total[[t.time]]*(rowSums(pm.a1[[t.time]]) -
                                            mu.a0.total[[t.time]]) + mu.a0.total[[t.time]]
    first_part <- w.cal[,1]*(y.imp - second_part1) # sum_part1 = sum(first_part)
    second_part <- second_part1 - mu.a0.total[[t.time]]

    # sum_part2 = mean(second_part)

    n.complete <- sum(r.mat[,t.time+1])
    index.complete <- which(r.mat[,t.time+1] == 1)

    pi.a0.obs <- matrix(0, nrow = n, ncol = t.time)
    pi.a1.obs <- matrix(0, nrow = n, ncol = t.time)
    ps.ratio.obs <- matrix(0, nrow = n, ncol = t.time)
    mu.a0.obs <- matrix(0, nrow = n, ncol = t.time+1)
    for(k in t.time:1){
      index <- which(r.mat[,k] == 1)
      pi.a0.obs[index,k] <- pi.a0.total[[t.time - k + 1]]
      pi.a1.obs[index,k] <- pi.a1.total[[t.time - k + 1]]
      ps.ratio.obs[index,k] <- ps.ratio[[k]]
      mu.a0.obs[index,k] <- mu.a0.total[[t.time - k + 1]]
    }
    pi.a0.cum <- cbind(1, t(apply(pi.a0.obs, 1, cumprod)))
    pi.a0.cum <- ifelse(pi.a0.cum == 0, 2, pi.a0.cum)
    mu.a0.obs[index.complete, t.time + 1] <- data[index.complete,y.name[t.time]]
    pi.a1.obs <- cbind(1, pi.a1.obs)

    # recode: cumulative probability (using all observed),
    # sum.part2 / sum.part.seq
    sum.part2 <- rep(0, n)
    sum.part.seq <- list()
    sum.part.seq.long <- matrix(0, nrow = n, ncol = t.time)
    mu.diff <- matrix(0, nrow = n, ncol = t.time)
    for(k in 1:t.time){
      sum.part2 <- pi.a0.cum[,k]*(1 - pi.a1.obs[,k+1])*ps.ratio.obs[,k] + sum.part2
      sum.part.seq[[k]] <- sum.part2
      mu.diff[,k] <- mu.a0.obs[,t.time - k + 2] - mu.a0.obs[,t.time - k + 1]
    }

    sum.part4 <- rep(0, n)
    for(k in 1:t.time){
      fourth_part <- w.cal[,k+1]*(sum.part.seq[[k]] - 1)*mu.diff[,t.time - k + 1]
      sum.part4 <- fourth_part + sum.part4
    }
    fourth_part <- sum.part4

    point_value <- sum(first_part) + mean(second_part) + sum(fourth_part)
    var_value <- mean((n*first_part + second_part + n*fourth_part - point_value)^2)/n

    return(list(point_value = point_value,
                var_value = var_value))
  }

  est <- NULL
  ve <- NULL
  if("mr" %in% type){
    res <- mr_est(e.ps, pm.a1, mu.a0.total,
                  pi.a0.total, pi.a1.total, ps.ratio)
    est <- c(est, res$point_value)
    ve <- c(ve, res$var_value)
    names(ve)[length(ve)] <- "mr"
  }

  if("mr.norm" %in% type){
    res <- mr_alter_est(e.ps, pm.a1, mu.a0.total,
                        pi.a0.total, pi.a1.total, ps.ratio)
    est <- c(est, res$point_value)
    ve <- c(ve, res$var_value)
    names(ve)[length(ve)] <- "mr.norm"
  }

  if("mr.cal" %in% type){
    if(is.null(list.cal)){
      stop("For calibration-based estimators, need to specify `mat.cal`")
    }
    g_whole_mat <- data.matrix(list.cal[[1]])
    tilde_g <- colMeans(g_whole_mat)

    # (2) For ps
    g1_mat <- g_whole_mat[a == 1,]
    g0_mat <- g_whole_mat[a == 0,]
    index.a0 <- which(a == 0)

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

    w.cal <- matrix(0, nrow = n, ncol = t.time + 1)

    w1_cal <- weight_fn(g1_mat, tilde_g)
    w.cal[a == 1,1] <- w1_cal
    w0_cal <- weight_fn(g0_mat, tilde_g)
    w0.cal.total <- rep(0, n)
    w0.cal.total[index.a0] <- w0_cal
    wr.cal.marginal <- rep(1, n)
    for(k in 1:t.time){
      index.r <- which(r.mat[,k+1] == 1)
      gr_mat <- list.cal[[k+1]]
      tilde_g.adj <- colMeans(gr_mat[r.mat[,k] == 1,])
      wr.cal <- weight_fn(gr_mat[index.r,], tilde_g.adj)
      wr.cal.total <- rep(0, n)
      wr.cal.total[index.r] <- wr.cal
      wr.cal.marginal <- wr.cal.total*wr.cal.marginal
      w.cal[,k+1] <- w0.cal.total*wr.cal.marginal/sum(w0.cal.total*wr.cal.marginal)
    }
    res <- mr_weight_est(w.cal, pm.a1, mu.a0.total,pi.a0.total, pi.a1.total,
                         ps.ratio)
    est <- c(est, res$point_value)
    ve <- c(ve, res$var_value)
    names(ve)[length(ve)] <- "mr.cal"
  }

  if("rppm" %in% type){
    est <- c(est, rppm_est(pi.a1.total, pm.a1, mu.a0.total))
  }

  if("psom" %in% type){
    est <- c(est, psom_est(e.ps, mu.a0.total))
  }

  if("psom.norm" %in% type){
    est <- c(est, psom_alter_est(e.ps, mu.a0.total))
  }

  if("psrp" %in% type){
    est <- c(est, psrp_est(e.ps, pi.a0.total, pi.a1.total, ps.ratio))
  }

  if("psom.norm" %in% type){
    est <- c(est, psrp_alter_est(e.ps, pi.a0.total, pi.a1.total, ps.ratio))
  }
  names(est) <- type

  return(list(est = est,
              ve = ve))
}

#' @importFrom stats binomial glm lm na.omit predict rnorm
NULL
