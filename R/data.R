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


#' Simulated data from the longitudinal study
#'
#' The simulated dataset for the longitudinal study that includes two follow-up
#' visits.
#'
#' @format ## `dat2`
#' A data frame with 500 rows and 15 columns:
#' \describe{
#'   \item{x1, x2, x3, x4}{The original continuous baseline covariates that
#'   are used to generate the non-linear transformed covariates}
#'   \item{x5}{The original binary baseline covariate}
#'   \item{z1, z2, z3, z4}{The transformed continuous baseline covariates that
#'   are used to generate the data}
#'   \item{z5}{The transformed binary baseline covariate}
#'   \item{a}{The treatment assignment.}
#'   \item{r1}{The response status at the first visit, where `r1 = 1` indicates
#'   the subject is observed at the first visit, and `r1 = 0`` implies the
#'   subject is missing at the first visit.}
#'   \item{r2}{The response status at the second visit, where `r2 = 1` indicates
#'   the subject is observed at the second visit, and `r2 = 0` implies the
#'   subject is missing at the second visit.}
#'   \item{y1}{The outcome at the first visit. Note that "NA" indicates that
#'   the subject is missing at the first visit, thus has no outcomes.}
#'   \item{y2}{The outcome at the second visit. Note that "NA" indicates that
#'   the subject is missing at the second visit, thus has no outcomes.}
#' }
"dat2"

#' An antidepressant trial data
#'
#' The data targets an antidepressant trial conducted under the Auspices of the
#' Drug Information Association. It was prepared by Mallinckrodt et al.(2014).
#' To illustrate the usage of our proposed estimator, we delete the individual
#' with intermittent missingness and only focus on the individuals with a
#' monotone missingness pattern, and we also delete three individuals with the
#' unobserved investigation site numbers. The reorganized dataset contain 99
#' subjects in the control group and 97 subjects in the treatment group. After
#' data preprocessing, 39 participants in the control group and 30 participants
#' in the treatment group suffered from monotone missingness.
#'
#' @format ## `datHAMD`
#' A data frame with 196 rows and 14 columns:
#' \describe{
#'   \item{id}{The original continuous baseline covariates that
#'   are used to generate the non-linear transformed covariates.}
#'   \item{trt}{The treatment assignment.}
#'   \item{poolinv}{The categorical variable indicating the investigation sites.}
#'   \item{base}{The baseline variable.}
#'   \item{y1,y2,y3,y4,y5}{The relative change from baseline at the k-th
#'   follow-up visit, with `NA` indicates the missing outcome and k = 1,...,5.}
#'   \item{r1,r2,r3,r4,r5}{The response status at the k-th follow-up visit,
#'   where `rk = 1` indicates the subject is observed at the k-th visit, and
#'   `rk = 0` implies the subject is missing at the k-th visit, for k = 1,...,5.}
#' }
#' @source <https://www.lshtm.ac.uk/research/centres-projects-groups/missing-data#dia-missing-data>
"datHAMD"
