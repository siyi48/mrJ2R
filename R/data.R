#' Simulated data from the cross-sectional study
#'
#' The simulated dataset for the cross-sectional study, i.e., the study with one
#' follow-up visit.
#'
#' @format ## `dat1`
#' A data frame with 500 rows and 13 columns:
#' \describe{
#'   \item{x1, x2, x3, x4}{The original continuous baseline covariates that
#'   are used to generate the non-linear transformed covariates}
#'   \item{x5}{The original binary baseline covariate}
#'   \item{z1, z2, z3, z4}{The transformed continuous baseline covariates that
#'   are used to generate the data}
#'   \item{z5}{The transformed binary baseline covariate}
#'   \item{a}{The treatment assignment.}
#'   \item{r}{The response status, where "r = 1" indicates the subject completes
#'   the study, and "r = 0" implies the subject is missing}
#'   \item{y}{The outcome. Note that "NA" indicates that the subject is missing,
#'   thus has no outcomes.}
#' }
"dat1"