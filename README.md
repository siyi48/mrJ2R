
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Codes Description

<!-- badges: start -->
<!-- badges: end -->

The goal of mrJ2R is to obtain the proposed estimators to evaluate the
average treatment effect (ATE) under jump-to-reference (J2R) based on
the paper **Multiply robust estimators in longitudinal studies with
missing data under control-based imputation**.

The file `sim_1time.R` contains four functions to get estimators under
the cross-sectional study as follows.

- `sim_dat(n, k, alpha, gamma)`: generate the simulated dataset;
- `model_est(dat)`: estimate the ATE using eight proposed estimators,
  including the three triply robust estimators $\hat \tau_{\text{tr}}$,
  $\hat \tau_{\text{tr-N}}$, and $\hat \tau_{\text{tr-C}}$ and five
  simple estimators $\hat \tau_{\text{ps-rp}}$,
  $\hat \tau_{\text{ps-rp-N}}$, $\hat \tau_{\text{ps-om}}$,
  $\hat \tau_{\text{ps-om-N}}$, and $\hat \tau_{\text{rp-om}}$;
- `nonpara_fn(dat, B, psrpom_est, psom_est, ps_est, rpom_est, rp_est, psrp_est, om_est, none_est)`:
  estimate the variation of the estimators using nonparmetric bootstrap,
  symmetric-t bootstrap, and bootstrap percentile;
- `main(seed)`: return the point and variance estimation results.

The file `sim_longi.R` contains four functions to get estimators under
the longitudinal study as follows.

- `sim_dat(n, k, alpha, gamma1, gamma2)`: generate the simulated
  dataset;
- `model_est(dat)`: estimate the ATE using eight proposed estimators,
  including the three multiply robust estimators
  $\hat \tau_{\text{mr}}$, $\hat \tau_{\text{mr-N}}$, and
  $\hat \tau_{\text{mr-C}}$ and five simple estimators
  $\hat \tau_{\text{ps-rp}}$, $\hat \tau_{\text{ps-rp-N}}$,
  $\hat \tau_{\text{ps-om}}$, $\hat \tau_{\text{ps-om-N}}$, and
  $\hat \tau_{\text{rp-pm}}$;
- `nonpara_fn(dat, B, point_est)`: estimate the variation of the
  estimators using nonparmetric bootstrap, symmetric-t bootstrap, and
  bootstrap percentile;
- `main(seed)`: return the point and variance estimation results.
