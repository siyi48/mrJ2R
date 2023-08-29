## For the cross-sectional studies ----
figure_output <- function(res){
  nsim <- length(res)
  psrpom_all <- matrix(0, 8, nsim)
  psom_all <- matrix(0, 8, nsim)
  ps_all <- matrix(0, 8, nsim)
  rpom_all <- matrix(0, 8, nsim)
  rp_all <- matrix(0, 8, nsim)
  psrp_all <- matrix(0, 8, nsim)
  om_all <- matrix(0, 8, nsim)
  none_all <- matrix(0, 8, nsim)
  true_rpom <- c(0)
  for(s in 1:nsim){
    psrpom_all[,s] <- res[[s]]$psrpom_all
    psom_all[,s] <- res[[s]]$psom_all
    ps_all[,s] <- res[[s]]$ps_all
    rpom_all[,s] <- res[[s]]$rpom_all
    rp_all[,s] <- res[[s]]$rp_all
    psrp_all[,s] <- res[[s]]$psrp_all
    om_all[,s] <- res[[s]]$om_all
    none_all[,s] <- res[[s]]$none_all
    true_rpom[s] <- res[[s]]$true_value_rpom
  }
  true_value <- mean(true_rpom)
  mat_all <- cbind(t(psrpom_all),
                   t(psom_all),
                   t(ps_all),
                   t(rpom_all),
                   t(rp_all),
                   t(psrp_all),
                   t(om_all),
                   t(none_all)) - true_value
  name_mat <- c(rep("tr",nsim),
                rep("tr-N",nsim),
                # rep("TR-GAM",nsim),
                # rep("TR-GAM2",nsim),
                rep("tr-C",nsim),
                rep("psrp",nsim),
                rep("psrp-N",nsim),
                rep("psom",nsim),
                rep("psom-N",nsim),
                rep("rpom",nsim))
  name_all <- rep(name_mat, 8)
  name_all <- factor(name_mat, levels = c("tr", "tr-N", "tr-C",
                                          "psrp", "psrp-N",
                                          "psom", "psom-N",
                                          "rpom"))
  ps_model <- c(rep("ps: yes", 3*8*nsim),
                rep("ps: no", 2*8*nsim),
                rep("ps: yes", 8*nsim),
                rep("ps: no", 2*8*nsim))
  ps_model <- factor(ps_model, levels = c("ps: yes", "ps: no"))
  rp_model <- c(rep("rp: yes", 8*nsim),
                rep("rp: no", 2*8*nsim),
                rep("rp: yes", 3*8*nsim),
                rep("rp: no", 2*8*nsim))
  rp_model <- factor(rp_model, levels = c("rp: yes", "rp: no"))
  om_model <- c(rep("om: yes", 2*8*nsim),
                rep("om: no", 1*8*nsim),
                rep("om: yes", 8*nsim),
                rep("om: no", 2*8*nsim),
                rep("om: yes", 8*nsim),
                rep("om: no", 8*nsim))
  om_model <- factor(om_model, levels = c("om: yes", "om: no"))

  res_df_all <- data.frame(bias = as.vector(mat_all),
                           est_name = name_all,
                           ps_model = ps_model,
                           rp_model = rp_model,
                           om_model = om_model)
  trc_ind <- with(res_df_all, ifelse(est_name == "tr-C", 1, 0))
  res_df_all$trc_ind <- trc_ind
  My_Theme = theme(
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13),
    legend.position = "none")
  p <- ggplot(res_df_all, aes(x=est_name, y=bias, color = trc_ind)) +
    facet_grid(ps_model ~ rp_model + om_model) +
    geom_boxplot() + geom_hline(yintercept=0, color = "grey") +
    xlab("estimator name") +
    ylim(-0.8,0.5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  p + My_Theme
}

assess_function <- function(res){
  nsim <- length(res)
  psrpom_all <- matrix(0, 8, nsim)
  psom_all <- matrix(0, 8, nsim)
  ps_all <- matrix(0, 8, nsim)
  rpom_all <- matrix(0, 8, nsim)
  rp_all <- matrix(0, 8, nsim)
  psrp_all <- matrix(0, 8, nsim)
  om_all <- matrix(0, 8, nsim)
  none_all <- matrix(0, 8, nsim)
  true_rpom <- c(0)
  # var estimate (nonparametric bootstrap)
  var_psrpom_all_np <- matrix(0, 8, nsim)
  var_psom_all_np <- matrix(0, 8, nsim)
  var_ps_all_np <- matrix(0, 8, nsim)
  var_rpom_all_np <- matrix(0, 8, nsim)
  var_rp_all_np <- matrix(0, 8, nsim)
  var_psrp_all_np <- matrix(0, 8, nsim)
  var_om_all_np <- matrix(0, 8, nsim)
  var_none_all_np <- matrix(0, 8, nsim)
  # var estimate (asymptotic theory)
  var_psrpom_all_asym <- matrix(0, 3, nsim)
  var_psom_all_asym <- matrix(0, 3, nsim)
  var_ps_all_asym <- matrix(0, 3, nsim)
  var_rpom_all_asym <- matrix(0, 3, nsim)
  var_rp_all_asym <- matrix(0, 3, nsim)
  var_psrp_all_asym <- matrix(0, 3, nsim)
  var_om_all_asym <- matrix(0, 3, nsim)
  var_none_all_asym <- matrix(0, 3, nsim)
  # c-star
  c_psrpom_star <- matrix(0, 3, nsim)
  c_psom_star <- matrix(0, 3, nsim)
  c_ps_star <- matrix(0, 3, nsim)
  c_rpom_star <- matrix(0, 3, nsim)
  c_rp_star <- matrix(0, 3, nsim)
  c_psrp_star <- matrix(0, 3, nsim)
  c_om_star <- matrix(0, 3, nsim)
  c_none_star <- matrix(0, 3, nsim)
  for(s in 1:nsim){
    psrpom_all[,s] <- res[[s]]$psrpom_all
    psom_all[,s] <- res[[s]]$psom_all
    ps_all[,s] <- res[[s]]$ps_all
    rpom_all[,s] <- res[[s]]$rpom_all
    rp_all[,s] <- res[[s]]$rp_all
    psrp_all[,s] <- res[[s]]$psrp_all
    om_all[,s] <- res[[s]]$om_all
    none_all[,s] <- res[[s]]$none_all
    true_rpom[s] <- res[[s]]$true_value_rpom
    # var estimates (nonparametric bootstrap)
    var_psrpom_all_np[,s] <- res[[s]]$var_psrpom_all_np
    var_psom_all_np[,s] <- res[[s]]$var_psom_all_np
    var_ps_all_np[,s] <- res[[s]]$var_ps_all_np
    var_rpom_all_np[,s] <- res[[s]]$var_rpom_all_np
    var_rp_all_np[,s] <- res[[s]]$var_rp_all_np
    var_psrp_all_np[,s] <- res[[s]]$var_psrp_all_np
    var_om_all_np[,s] <- res[[s]]$var_om_all_np
    var_none_all_np[,s] <- res[[s]]$var_none_all_np
    # var estimates (asymptotic theory)
    var_psrpom_all_asym[,s] <- res[[s]]$var_psrpom_all_asym
    var_psom_all_asym[,s] <- res[[s]]$var_psom_all_asym
    var_ps_all_asym[,s] <- res[[s]]$var_ps_all_asym
    var_rpom_all_asym[,s] <- res[[s]]$var_rpom_all_asym
    var_rp_all_asym[,s] <- res[[s]]$var_rp_all_asym
    var_psrp_all_asym[,s] <- res[[s]]$var_psrp_all_asym
    var_om_all_asym[,s] <- res[[s]]$var_om_all_asym
    var_none_all_asym[,s] <- res[[s]]$var_none_all_asym
    # c-star (symmetric t)
    c_psrpom_star[,s] <- res[[s]]$c_psrpom_star
    c_psom_star[,s] <- res[[s]]$c_psom_star
    c_ps_star[,s] <- res[[s]]$c_ps_star
    c_rpom_star[,s] <- res[[s]]$c_rpom_star
    c_rp_star[,s] <- res[[s]]$c_rp_star
    c_psrp_star[,s] <- res[[s]]$c_psrp_star
    c_om_star[,s] <- res[[s]]$c_om_star
    c_none_star[,s] <- res[[s]]$c_none_star
  }
  true_value <- mean(true_rpom)

  mean_psrpom_all <- apply(psrpom_all, 1, mean)
  truevar_psrpom_all <- apply(psrpom_all, 1, var)
  bootvar_psrpom_all <- apply(var_psrpom_all_np, 1, mean)
  asymvar_psrpom_all <- apply(var_psrpom_all_asym, 1, mean)

  mean_psom_all <- apply(psom_all, 1, mean)
  truevar_psom_all <- apply(psom_all, 1, var)
  bootvar_psom_all <- apply(var_psom_all_np, 1, mean)
  asymvar_psom_all <- apply(var_psom_all_asym, 1, mean)

  mean_ps_all <- apply(ps_all, 1, mean)
  truevar_ps_all <- apply(ps_all, 1, var)
  bootvar_ps_all <- apply(var_ps_all_np, 1, mean)
  asymvar_ps_all <- apply(var_ps_all_asym, 1, mean)

  mean_rpom_all <- apply(rpom_all, 1, mean)
  truevar_rpom_all <- apply(rpom_all, 1, var)
  bootvar_rpom_all <- apply(var_rpom_all_np, 1, mean)
  asymvar_rpom_all <- apply(var_rpom_all_asym, 1, mean)

  mean_rp_all <- apply(rp_all, 1, mean)
  truevar_rp_all <- apply(rp_all, 1, var)
  bootvar_rp_all <- apply(var_rp_all_np, 1, mean)
  asymvar_rp_all <- apply(var_rp_all_asym, 1, mean)

  mean_psrp_all <- apply(psrp_all, 1, mean)
  truevar_psrp_all <- apply(psrp_all, 1, var)
  bootvar_psrp_all <- apply(var_psrp_all_np, 1, mean)
  asymvar_psrp_all <- apply(var_psrp_all_asym, 1, mean)

  mean_om_all <- apply(om_all, 1, mean)
  truevar_om_all <- apply(om_all, 1, var)
  bootvar_om_all <- apply(var_om_all_np, 1, mean)
  asymvar_om_all <- apply(var_om_all_asym, 1, mean)

  mean_none_all <- apply(none_all, 1, mean)
  truevar_none_all <- apply(none_all, 1, var)
  bootvar_none_all <- apply(var_none_all_np, 1, mean)
  asymvar_none_all <- apply(var_none_all_asym, 1, mean)

  # relative bias
  rela_psrpom <- (bootvar_psrpom_all - truevar_psrpom_all)/truevar_psrpom_all
  rela_psom <- (bootvar_psom_all - truevar_psom_all)/truevar_psom_all
  rela_ps <- (bootvar_ps_all - truevar_ps_all)/truevar_ps_all
  rela_rpom <- (bootvar_rpom_all - truevar_rpom_all)/truevar_rpom_all
  rela_rp <- (bootvar_rp_all - truevar_rp_all)/truevar_rp_all
  rela_psrp <- (bootvar_psrp_all - truevar_psrp_all)/truevar_psrp_all
  rela_om <- (bootvar_om_all - truevar_om_all)/truevar_om_all
  rela_none <- (bootvar_none_all - truevar_none_all)/truevar_none_all

  # coverage rate
  ## Case 1
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- psrpom_all - 1.96*sqrt(var_psrpom_all_np)
  ci_upper <- psrpom_all + 1.96*sqrt(var_psrpom_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psrpom_np <- apply(cover_mat, 1, mean)
  ci_width_psrpom_np <- rowMeans(ci_upper - ci_lower)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- psrpom_all[1:3,] - 1.96*sqrt(var_psrpom_all_asym)
  ci_upper <- psrpom_all[1:3,] + 1.96*sqrt(var_psrpom_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psrpom_asym <- apply(cover_mat, 1, mean)
  ci_width_psrpom_asym <- rowMeans(ci_upper - ci_lower)
  ## (3) Symmetric t
  ci_lower <- psrpom_all[1:3,] - c_psrpom_star*sqrt(var_psrpom_all_asym)
  ci_upper <- psrpom_all[1:3,] + c_psrpom_star*sqrt(var_psrpom_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psrpom_symmt <- apply(cover_mat, 1, mean)
  ci_width_psrpom_symmt <- rowMeans(ci_upper - ci_lower)

  ## Case 2
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- psom_all - 1.96*sqrt(var_psom_all_np)
  ci_upper <- psom_all + 1.96*sqrt(var_psom_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psom_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- psom_all[1:3,] - 1.96*sqrt(var_psom_all_asym)
  ci_upper <- psom_all[1:3,] + 1.96*sqrt(var_psom_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psom_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- psom_all[1:3,] - c_psom_star*sqrt(var_psom_all_asym)
  ci_upper <- psom_all[1:3,] + c_psom_star*sqrt(var_psom_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psom_symmt <- apply(cover_mat, 1, mean)

  ## Case 3
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- ps_all - 1.96*sqrt(var_ps_all_np)
  ci_upper <- ps_all + 1.96*sqrt(var_ps_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_ps_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- ps_all[1:3,] - 1.96*sqrt(var_ps_all_asym)
  ci_upper <- ps_all[1:3,] + 1.96*sqrt(var_ps_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_ps_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- ps_all[1:3,] - c_ps_star*sqrt(var_ps_all_asym)
  ci_upper <- ps_all[1:3,] + c_ps_star*sqrt(var_ps_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_ps_symmt <- apply(cover_mat, 1, mean)

  ## Case 4
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- rpom_all - 1.96*sqrt(var_rpom_all_np)
  ci_upper <- rpom_all + 1.96*sqrt(var_rpom_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_rpom_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- rpom_all[1:3,] - 1.96*sqrt(var_rpom_all_asym)
  ci_upper <- rpom_all[1:3,] + 1.96*sqrt(var_rpom_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_rpom_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- rpom_all[1:3,] - c_rpom_star*sqrt(var_rpom_all_asym)
  ci_upper <- rpom_all[1:3,] + c_rpom_star*sqrt(var_rpom_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_rpom_symmt <- apply(cover_mat, 1, mean)

  ## Case 5
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- rp_all - 1.96*sqrt(var_rp_all_np)
  ci_upper <- rp_all + 1.96*sqrt(var_rp_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_rp_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- rp_all[1:3,] - 1.96*sqrt(var_rp_all_asym)
  ci_upper <- rp_all[1:3,] + 1.96*sqrt(var_rp_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_rp_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- rp_all[1:3,] - c_rp_star*sqrt(var_rp_all_asym)
  ci_upper <- rp_all[1:3,] + c_rp_star*sqrt(var_rp_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_rp_symmt <- apply(cover_mat, 1, mean)

  ## Case 6
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- psrp_all - 1.96*sqrt(var_psrp_all_np)
  ci_upper <- psrp_all + 1.96*sqrt(var_psrp_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psrp_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- psrp_all[1:3,] - 1.96*sqrt(var_psrp_all_asym)
  ci_upper <- psrp_all[1:3,] + 1.96*sqrt(var_psrp_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psrp_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- psrp_all[1:3,] - c_psrp_star*sqrt(var_psrp_all_asym)
  ci_upper <- psrp_all[1:3,] + c_psrp_star*sqrt(var_psrp_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_psrp_symmt <- apply(cover_mat, 1, mean)

  ## Case 7
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- om_all - 1.96*sqrt(var_om_all_np)
  ci_upper <- om_all + 1.96*sqrt(var_om_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_om_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- om_all[1:3,] - 1.96*sqrt(var_om_all_asym)
  ci_upper <- om_all[1:3,] + 1.96*sqrt(var_om_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_om_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- om_all[1:3,] - c_om_star*sqrt(var_om_all_asym)
  ci_upper <- om_all[1:3,] + c_om_star*sqrt(var_om_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_om_symmt <- apply(cover_mat, 1, mean)

  ## Case 8
  ## (1) Wald-type (nonparametric bootstrap)
  ci_lower <- none_all - 1.96*sqrt(var_none_all_np)
  ci_upper <- none_all + 1.96*sqrt(var_none_all_np)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_none_np <- apply(cover_mat, 1, mean)
  ## (2) Wald-type (asymptotic variance)
  ci_lower <- none_all[1:3,] - 1.96*sqrt(var_none_all_asym)
  ci_upper <- none_all[1:3,] + 1.96*sqrt(var_none_all_asym)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_none_asym <- apply(cover_mat, 1, mean)
  ## (3) Symmetric t
  ci_lower <- none_all[1:3,] - c_none_star*sqrt(var_none_all_asym)
  ci_upper <- none_all[1:3,] + c_none_star*sqrt(var_none_all_asym)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_none_symmt <- apply(cover_mat, 1, mean)

  # bootstrap percentile
  cover_rate_psrpom_bprec <- rep(0,8)
  ci_width_psrpom_bprec <- rep(0,8)
  cover_rate_psom_bprec <- rep(0,8)
  ci_width_psom_bprec <- rep(0,8)
  cover_rate_ps_bprec <- rep(0,8)
  ci_width_ps_bprec <- rep(0,8)
  cover_rate_rpom_bprec <- rep(0,8)
  ci_width_rpom_bprec <- rep(0,8)
  cover_rate_rp_bprec <- rep(0,8)
  ci_width_rp_bprec <- rep(0,8)
  cover_rate_psrp_bprec <- rep(0,8)
  ci_width_psrp_bprec <- rep(0,8)
  cover_rate_om_bprec <- rep(0,8)
  ci_width_om_bprec <- rep(0,8)
  cover_rate_none_bprec <- rep(0,8)
  ci_width_none_bprec <- rep(0,8)
  for(s in 1:nsim){
    ## Case 1
    percentile_psrpom_boot <- res[[s]]$percentile_psrpom_boot
    cover_rate_psrpom_bprec <- apply(percentile_psrpom_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_psrpom_bprec
    ci_width_psrpom_bprec <- apply(percentile_psrpom_boot, 2, function(x)
      x[2] - x[1]) + ci_width_psrpom_bprec
    ## Case 2
    percentile_psom_boot <- res[[s]]$percentile_psom_boot
    cover_rate_psom_bprec <- apply(percentile_psom_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_psom_bprec
    ci_width_psom_bprec <- apply(percentile_psom_boot, 2, function(x)
      x[2] - x[1]) + ci_width_psom_bprec
    ## Case 3
    percentile_ps_boot <- res[[s]]$percentile_ps_boot
    cover_rate_ps_bprec <- apply(percentile_ps_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_ps_bprec
    ci_width_ps_bprec <- apply(percentile_ps_boot, 2, function(x)
      x[2] - x[1]) + ci_width_ps_bprec
    ## Case 4
    percentile_rpom_boot <- res[[s]]$percentile_rpom_boot
    cover_rate_rpom_bprec <- apply(percentile_rpom_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_rpom_bprec
    ci_width_rpom_bprec <- apply(percentile_rpom_boot, 2, function(x)
      x[2] - x[1]) + ci_width_rpom_bprec
    ## Case 5
    percentile_rp_boot <- res[[s]]$percentile_rp_boot
    cover_rate_rp_bprec <- apply(percentile_rp_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_rp_bprec
    ci_width_rp_bprec <- apply(percentile_rp_boot, 2, function(x)
      x[2] - x[1]) + ci_width_rp_bprec
    ## Case 6
    percentile_psrp_boot <- res[[s]]$percentile_psrp_boot
    cover_rate_psrp_bprec <- apply(percentile_psrp_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_psrp_bprec
    ci_width_psrp_bprec <- apply(percentile_psrp_boot, 2, function(x)
      x[2] - x[1]) + ci_width_psrp_bprec
    ## Case 7
    percentile_om_boot <- res[[s]]$percentile_om_boot
    cover_rate_om_bprec <- apply(percentile_om_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_om_bprec
    ci_width_om_bprec <- apply(percentile_om_boot, 2, function(x)
      x[2] - x[1]) + ci_width_om_bprec
    ## Case 8
    percentile_none_boot <- res[[s]]$percentile_none_boot
    cover_rate_none_bprec <- apply(percentile_none_boot, 2, function(x)
      ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_none_bprec
    ci_width_none_bprec <- apply(percentile_none_boot, 2, function(x)
      x[2] - x[1]) + ci_width_none_bprec
  }
  cover_rate_psrpom_bprec <- cover_rate_psrpom_bprec/nsim
  ci_width_psrpom_bprec <- ci_width_psrpom_bprec/nsim
  cover_rate_psom_bprec <- cover_rate_psom_bprec/nsim
  ci_width_psom_bprec <- ci_width_psom_bprec/nsim
  cover_rate_ps_bprec <- cover_rate_ps_bprec/nsim
  ci_width_ps_bprec <- ci_width_ps_bprec/nsim
  cover_rate_rpom_bprec <- cover_rate_rpom_bprec/nsim
  ci_width_rpom_bprec <- ci_width_rpom_bprec/nsim
  cover_rate_rp_bprec <- cover_rate_rp_bprec/nsim
  ci_width_rp_bprec <- ci_width_rp_bprec/nsim
  cover_rate_psrp_bprec <- cover_rate_psrp_bprec/nsim
  ci_width_psrp_bprec <- ci_width_psrp_bprec/nsim
  cover_rate_om_bprec <- cover_rate_om_bprec/nsim
  ci_width_om_bprec <- ci_width_om_bprec/nsim
  cover_rate_none_bprec <- cover_rate_none_bprec/nsim
  ci_width_none_bprec <- ci_width_none_bprec/nsim

  ## Table (only for psrpom + 3 TR etsimators)
  table_mat <- cbind(cover_rate_psrpom_np[1:3],
                     ci_width_psrpom_np[1:3],
                     cover_rate_psrpom_asym,
                     ci_width_psrpom_asym,
                     cover_rate_psrpom_symmt,
                     ci_width_psrpom_symmt)*10^2

  return(list(# true value
    true_value = true_value,
    # psrpom
    mean_psrpom = mean_psrpom_all,
    truevar_psrpom = truevar_psrpom_all,
    bootvar_psrpom = bootvar_psrpom_all,
    asymvar_psrpom = asymvar_psrpom_all,
    rela_psrpom = rela_psrpom,
    cover_rate_psrpom_np = cover_rate_psrpom_np,
    cover_rate_psrpom_asym = cover_rate_psrpom_asym,
    cover_rate_psrpom_symmt = cover_rate_psrpom_symmt,
    cover_rate_psrpom_bprec = cover_rate_psrpom_bprec,
    # psom
    mean_psom = mean_psom_all,
    truevar_psom = truevar_psom_all,
    bootvar_psom = bootvar_psom_all,
    asymvar_psom = asymvar_psom_all,
    rela_psom = rela_psom,
    cover_rate_psom_np = cover_rate_psom_np,
    cover_rate_psom_asym = cover_rate_psom_asym,
    cover_rate_psom_symmt = cover_rate_psom_symmt,
    cover_rate_psom_bprec = cover_rate_psom_bprec,
    # ps
    mean_ps = mean_ps_all,
    truevar_ps = truevar_ps_all,
    bootvar_ps = bootvar_ps_all,
    asymvar_ps = asymvar_ps_all,
    rela_ps = rela_ps,
    cover_rate_ps_np = cover_rate_ps_np,
    cover_rate_ps_asym = cover_rate_ps_asym,
    cover_rate_ps_symmt = cover_rate_ps_symmt,
    cover_rate_ps_bprec = cover_rate_ps_bprec,
    # rpom
    mean_rpom = mean_rpom_all,
    truevar_rpom = truevar_rpom_all,
    bootvar_rpom = bootvar_rpom_all,
    asymvar_rpom = asymvar_rpom_all,
    rela_rpom = rela_rpom,
    cover_rate_rpom_np = cover_rate_rpom_np,
    cover_rate_rpom_asym = cover_rate_rpom_asym,
    cover_rate_rpom_symmt = cover_rate_rpom_symmt,
    cover_rate_rpom_bprec = cover_rate_rpom_bprec,
    # rp
    mean_rp = mean_rp_all,
    truevar_rp = truevar_rp_all,
    bootvar_rp = bootvar_rp_all,
    asymvar_rp = asymvar_rp_all,
    rela_rp = rela_rp,
    cover_rate_rp_np = cover_rate_rp_np,
    cover_rate_rp_asym = cover_rate_rp_asym,
    cover_rate_rp_symmt = cover_rate_rp_symmt,
    cover_rate_rp_bprec = cover_rate_rp_bprec,
    # psrp
    mean_psrp = mean_psrp_all,
    truevar_psrp = truevar_psrp_all,
    bootvar_psrp = bootvar_psrp_all,
    asymvar_psrp = asymvar_psrp_all,
    rela_psrp = rela_psrp,
    cover_rate_psrp_np = cover_rate_psrp_np,
    cover_rate_psrp_asym = cover_rate_psrp_asym,
    cover_rate_psrp_symmt = cover_rate_psrp_symmt,
    cover_rate_psrp_bprec = cover_rate_psrp_bprec,
    # om
    mean_om = mean_om_all,
    truevar_om = truevar_om_all,
    bootvar_om = bootvar_om_all,
    asymvar_om = asymvar_om_all,
    rela_om = rela_om,
    cover_rate_om_np = cover_rate_om_np,
    cover_rate_om_asym = cover_rate_om_asym,
    cover_rate_om_symmt = cover_rate_om_symmt,
    cover_rate_om_bprec = cover_rate_om_bprec,
    # none
    mean_none = mean_none_all,
    truevar_none = truevar_none_all,
    bootvar_none = bootvar_none_all,
    asymvar_none = asymvar_none_all,
    rela_none = rela_none,
    cover_rate_none_np = cover_rate_none_np,
    cover_rate_none_asym = cover_rate_none_asym,
    cover_rate_none_symmt = cover_rate_none_symmt,
    cover_rate_none_bprec = cover_rate_none_bprec,
    # table
    table_mat = table_mat
  ))
}

table_function <- function(res){
  nsim <- length(res)
  psrpom_all <- matrix(0, 8, nsim)
  psom_all <- matrix(0, 8, nsim)
  ps_all <- matrix(0, 8, nsim)
  rpom_all <- matrix(0, 8, nsim)
  rp_all <- matrix(0, 8, nsim)
  psrp_all <- matrix(0, 8, nsim)
  om_all <- matrix(0, 8, nsim)
  none_all <- matrix(0, 8, nsim)
  true_rpom <- c(0)
  var_psrpom_all_np <- matrix(0, 8, nsim)
  var_psom_all_np <- matrix(0, 8, nsim)
  var_ps_all_np <- matrix(0, 8, nsim)
  var_rpom_all_np <- matrix(0, 8, nsim)
  var_rp_all_np <- matrix(0, 8, nsim)
  var_psrp_all_np <- matrix(0, 8, nsim)
  var_om_all_np <- matrix(0, 8, nsim)
  var_none_all_np <- matrix(0, 8, nsim)
  for(s in 1:nsim){
    psrpom_all[,s] <- res[[s]]$psrpom_all
    psom_all[,s] <- res[[s]]$psom_all
    ps_all[,s] <- res[[s]]$ps_all
    rpom_all[,s] <- res[[s]]$rpom_all
    rp_all[,s] <- res[[s]]$rp_all
    psrp_all[,s] <- res[[s]]$psrp_all
    om_all[,s] <- res[[s]]$om_all
    none_all[,s] <- res[[s]]$none_all
    true_rpom[s] <- res[[s]]$true_value_rpom
    var_psrpom_all_np[,s] <- res[[s]]$var_psrpom_all_np
    var_psom_all_np[,s] <- res[[s]]$var_psom_all_np
    var_ps_all_np[,s] <- res[[s]]$var_ps_all_np
    var_rpom_all_np[,s] <- res[[s]]$var_rpom_all_np
    var_rp_all_np[,s] <- res[[s]]$var_rp_all_np
    var_psrp_all_np[,s] <- res[[s]]$var_psrp_all_np
    var_om_all_np[,s] <- res[[s]]$var_om_all_np
    var_none_all_np[,s] <- res[[s]]$var_none_all_np
  }
  true_value <- mean(true_rpom)

  mean_psrpom_all <- apply(psrpom_all, 1, mean)
  truevar_psrpom_all <- apply(psrpom_all, 1, var)
  bootvar_psrpom_all <- apply(var_psrpom_all_np, 1, mean)

  round((mean_psrpom_all - true_value)*100, digits = 2)
  round(sqrt(truevar_psrpom_all)*100, digits = 2)

  ci_mat <- cbind(mean_psrpom_all-1.96*sqrt(bootvar_psrpom_all),
                  mean_psrpom_all+1.96*sqrt(bootvar_psrpom_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)


  mean_psom_all <- apply(psom_all, 1, mean)
  truevar_psom_all <- apply(psom_all, 1, var)
  bootvar_psom_all <- apply(var_psom_all_np, 1, mean)
  round((mean_psom_all - true_value)*100, digits = 2)
  round(sqrt(truevar_psom_all)*100, digits = 2)

  ci_mat <- cbind(mean_psom_all-1.96*sqrt(bootvar_psom_all),
                  mean_psom_all+1.96*sqrt(bootvar_psom_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)

  mean_ps_all <- apply(ps_all, 1, mean)
  truevar_ps_all <- apply(ps_all, 1, var)
  bootvar_ps_all <- apply(var_ps_all_np, 1, mean)
  round((mean_ps_all - true_value)*100, digits = 2)
  round(sqrt(truevar_ps_all)*100, digits = 2)

  ci_mat <- cbind(mean_ps_all-1.96*sqrt(bootvar_ps_all),
                  mean_ps_all+1.96*sqrt(bootvar_ps_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)

  mean_rpom_all <- apply(rpom_all, 1, mean)
  truevar_rpom_all <- apply(rpom_all, 1, var)
  bootvar_rpom_all <- apply(var_rpom_all_np, 1, mean)
  round((mean_rpom_all - true_value)*100, digits = 2)
  round(sqrt(truevar_rpom_all)*100, digits = 2)

  ci_mat <- cbind(mean_rpom_all-1.96*sqrt(bootvar_rpom_all),
                  mean_rpom_all+1.96*sqrt(bootvar_rpom_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)

  mean_rp_all <- apply(rp_all, 1, mean)
  truevar_rp_all <- apply(rp_all, 1, var)
  bootvar_rp_all <- apply(var_rp_all_np, 1, mean)
  round((mean_rp_all - true_value)*100, digits = 2)
  round(sqrt(truevar_rp_all)*100, digits = 2)

  ci_mat <- cbind(mean_rp_all-1.96*sqrt(bootvar_rp_all),
                  mean_rp_all+1.96*sqrt(bootvar_rp_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)

  mean_psrp_all <- apply(psrp_all, 1, mean)
  truevar_psrp_all <- apply(psrp_all, 1, var)
  bootvar_psrp_all <- apply(var_psrp_all_np, 1, mean)
  round((mean_psrp_all - true_value)*100, digits = 2)
  round(sqrt(truevar_psrp_all)*100, digits = 2)


  ci_mat <- cbind(mean_psrp_all-1.96*sqrt(bootvar_psrp_all),
                  mean_psrp_all+1.96*sqrt(bootvar_psrp_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)

  mean_om_all <- apply(om_all, 1, mean)
  truevar_om_all <- apply(om_all, 1, var)
  bootvar_om_all <- apply(var_om_all_np, 1, mean)
  round((mean_om_all - true_value)*100, digits = 2)
  round(sqrt(truevar_om_all)*100, digits = 2)

  ci_mat <- cbind(mean_om_all-1.96*sqrt(bootvar_om_all),
                  mean_om_all+1.96*sqrt(bootvar_om_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)

  mean_none_all <- apply(none_all, 1, mean)
  truevar_none_all <- apply(none_all, 1, var)
  bootvar_none_all <- apply(var_none_all_np, 1, mean)
  round((mean_none_all - true_value)*100, digits = 2)
  round(sqrt(truevar_none_all)*100, digits = 2)

  ci_mat <- cbind(mean_none_all-1.96*sqrt(bootvar_none_all),
                  mean_none_all+1.96*sqrt(bootvar_none_all))
  ci_length <- ci_mat[,2] - ci_mat[,1]
  round(ci_length*100, digits = 2)
}

## For the longitudinal studies ----
figure_output <- function(res){
  nsim <- length(res)
  psrpom_all <- matrix(0, 8, nsim)
  true_psom <- c(0)
  for(s in 1:nsim){
    psrpom_all[,s] <- res[[s]]$point_est
    true_psom[s] <- res[[s]]$true_value
  }
  true_value <- mean(true_psom)
  mat_all <- t(psrpom_all) - true_value
  name_mat <- c(rep("mr",nsim),
                rep("mr-N",nsim),
                rep("mr-C",nsim),
                rep("psrp",nsim),
                rep("psrp-N",nsim),
                rep("psom",nsim),
                rep("psom-N",nsim),
                rep("rppm",nsim))
  name_all <- rep(name_mat, 8)
  name_all <- factor(name_mat, levels = c("mr", "mr-N", "mr-C",
                                          "psrp", "psrp-N",
                                          "psom", "psom-N",
                                          "rppm"))

  res_df_all <- data.frame(bias = as.vector(mat_all),
                           est_name = name_all)
  trc_ind <- with(res_df_all, ifelse(est_name == "mr-C", 1, 0))
  res_df_all$trc_ind <- trc_ind
  My_Theme = theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    legend.position = "none")
  p <- ggplot(res_df_all, aes(x=est_name, y=bias, color = trc_ind)) +
    geom_boxplot() + geom_hline(yintercept=0, color = "grey") +
    xlab("estimator name")
  # p <- ggplot(res_df_all[-which(abs(res_df_all$bias) >= 9),], aes(x=est_name, y=bias)) +
  #   facet_grid(ps_model ~ rp_model + om_model) +
  #   geom_boxplot() + geom_hline(yintercept=0, color = "grey") +
  #   xlab("estimator name")
  p + My_Theme
}

assess_function <- function(res){
  nsim <- length(res)
  psrpom_all <- matrix(0, 8, nsim)
  true_psom <- c(0)
  var_psrpom_all <- matrix(0, 8, nsim)
  c_star <- matrix(0, 8, nsim)
  var_boot <- matrix(0, 8, nsim)
  for(s in 1:nsim){
    psrpom_all[,s] <- res[[s]]$point_est
    true_psom[s] <- res[[s]]$true_value
    var_psrpom_all[,s] <- res[[s]]$var_est
    c_star[,s] <- res[[s]]$c_star
    var_boot[,s] <- res[[s]]$var_boot
  }
  true_value <- mean(true_psom)

  var_est_all <- rbind(var_psrpom_all[1:3,], var_boot[4:8,])

  mean_psrpom_all <- apply(psrpom_all, 1, mean)
  truevar_psrpom_all <- apply(psrpom_all, 1, var)
  bootvar_psrpom_all <- apply(var_est_all, 1, mean)

  # relative bias
  rela_psrpom <- (bootvar_psrpom_all - truevar_psrpom_all)/truevar_psrpom_all

  # coverage rate
  ci_lower <- psrpom_all - c_star*sqrt(var_est_all)
  ci_upper <- psrpom_all + c_star*sqrt(var_est_all)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_symmt <- apply(cover_mat, 1, mean)
  ci_width_symmt <- rowMeans(ci_upper - ci_lower)

  ci_lower <- psrpom_all - 1.96*sqrt(var_boot)
  ci_upper <- psrpom_all + 1.96*sqrt(var_boot)
  mat <- cbind(rowMeans(ci_lower), rowMeans(ci_upper))
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_np <- apply(cover_mat, 1, mean)
  ci_width_np <- rowMeans(ci_upper - ci_lower)

  cover_rate_bprec <- rep(0,8)
  ci_width_bprec <- rep(0,8)
  for(s in 1:nsim){
    percentile_boot <- res[[s]]$percentile_boot
    cover_rate_bprec <- apply(percentile_boot, 2, function(x) ifelse(x[1] < true_value & x[2] > true_value, 1, 0)) + cover_rate_bprec
    ci_width_bprec <- apply(percentile_boot, 2, function(x) x[2] - x[1]) + ci_width_bprec
  }
  cover_rate_bprec <- cover_rate_bprec/nsim
  ci_width_bprec <- ci_width_bprec/nsim

  ci_lower <- psrpom_all - 1.96*sqrt(var_est_all)
  ci_upper <- psrpom_all + 1.96*sqrt(var_est_all)
  cover_mat <- apply(ci_lower, 2, function(x) ifelse(x < true_value, 1, 0)) * apply(ci_upper, 2, function(x) ifelse(x > true_value, 1, 0))
  cover_rate_asym <- apply(cover_mat, 1, mean)
  ci_width_asym <- rowMeans(ci_upper - ci_lower)

  ## Table (only for 3 MR etsimators)
  table_mat_ci <- cbind(cover_rate_np[1:3],
                        ci_width_np[1:3],
                        cover_rate_asym[1:3],
                        ci_width_asym[1:3],
                        cover_rate_symmt[1:3],
                        ci_width_symmt[1:3])*10^2

  ## Table (all estimator)
  table_mat_all <- cbind(cover_rate_symmt, ci_width_symmt)*10^2

  return(list(# true value
    true_value = true_value,
    # psrpom
    mean_psrpom = mean_psrpom_all,
    truevar_psrpom = truevar_psrpom_all,
    bootvar_psrpom = bootvar_psrpom_all,
    rela_psrpom = rela_psrpom,
    # coverage
    cover_rate_symmt = cover_rate_symmt,
    cover_rate_bprec = cover_rate_bprec,
    cover_rate_asym = cover_rate_asym,
    cover_rate_np = cover_rate_np,
    # CI length
    ci_width_np = ci_width_np,
    ci_width_symmt = ci_width_symmt,
    ci_width_bprec = ci_width_bprec,
    ci_width_asym = ci_width_asym,
    # table
    table_mat_all = table_mat_all,
    table_mat_ci = table_mat_ci
  ))
}
