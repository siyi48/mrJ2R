#' Obtain the proposed estimator in the cross-sectional study
#'
#' This function implements the proposed RPSFT + Cox method and returns an
#' estimate of the treatment effect on the overall survival. It can be used if
#' investigators have sufficient knowledge of the PDT effects and the survival
#' time during the PDT period.
#'
#' @param formula.ps the matrix of all covariates included in the propensity
#' score model, with the dimension: n*p1, n: the sample size, p1: dimension of the
#' covariates.
#' @param formula.rp the matrix of all covariates included in the response status
#' model, with the dimension: n*p2, n: the sample size, p2: dimension of the
#' covariates.
#' @param formula.om the matrix of all covariates included in the outcome mean
#' model, with the dimension: n*p3, n: the sample size, p3: dimension of the
#' covariates.
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
#' \item `tr`: the triply robust estimator with calibration
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

#' @importFrom stats binomial glm lm na.omit predict rnorm
NULL
