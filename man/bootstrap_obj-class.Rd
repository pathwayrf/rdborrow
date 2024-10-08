% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_class.R
\docType{class}
\name{bootstrap_obj-class}
\alias{bootstrap_obj-class}
\alias{.bootstrap_obj}
\title{Bootstrap class}
\description{
Bootstrap class
}
\section{Slots}{

\describe{
\item{\code{replicates}}{}

\item{\code{bootstrap_CI_type}}{}
}}

\examples{
bootstrap_obj = setup_bootstrap(
  replicates = 2e3,
  bootstrap_CI_type = "perc"
)

method_weighting_obj = setup_method_weighting(
  method_name = "AIPW",
  optimal_weight_flag = T, 
  wt = 0,
  bootstrap_flag = T,
  bootstrap_obj = bootstrap_obj,
  model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5",
  model_form_mu0_ext = c("y1 ~ x1 + x2 + x3 + x4 + x5",
                         "y2 ~ x1 + x2 + x3 + x4 + x5"))

analysis_primary_obj = setup_analysis_primary(
  data = SyntheticData,
  trial_status = "S", 
  treatment = "A", 
  outcome = c("y1", "y2"), 
  covariates = c("x1", "x2", "x3", "x4", "x5"), 
  method_weighting_obj = method_weighting_obj)


res = run_analysis(analysis_primary_obj)
}
