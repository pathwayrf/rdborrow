% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DID_EC_AIPW_bootstrap.R
\name{DID_EC_AIPW_bootstrap}
\alias{DID_EC_AIPW_bootstrap}
\title{Difference in difference + AIPW + external control borrowing}
\usage{
DID_EC_AIPW_bootstrap(
  data,
  indices,
  outcome_col_name,
  trial_status_col_name,
  treatment_col_name,
  covariates_col_name,
  T_cross,
  model_form_piS = "",
  model_form_piA = "",
  model_form_mu0_ext = ""
)
}
\arguments{
\item{data}{}

\item{indices}{}

\item{outcome_col_name}{}

\item{trial_status_col_name}{}

\item{treatment_col_name}{}

\item{covariates_col_name}{}

\item{T_cross}{}

\item{model_form_piS}{}

\item{model_form_piA}{}

\item{model_form_mu0_ext}{}
}
\value{
tau and standard deviation
}
\description{
Difference in difference + AIPW + external control borrowing
}
\examples{
DID_EC_AIPW_bootstrap = function(data,
  indices,
  outcome_col_name, 
  trial_status_col_name, 
  treatment_col_name, 
  covariates_col_name,
  T_cross,
  model_form_piS = "",
  model_form_piA = "",
  model_form_mu0_ext = ""
)
}
