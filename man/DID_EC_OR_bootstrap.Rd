% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DID_EC_OR_bootstrap.R
\name{DID_EC_OR_bootstrap}
\alias{DID_EC_OR_bootstrap}
\title{Difference in difference + outcome regression + external control borrowing}
\usage{
DID_EC_OR_bootstrap(
  data = data,
  indices = indices,
  outcome_col_name,
  trial_status_col_name,
  treatment_col_name,
  covariates_col_name,
  T_cross,
  model_form_mu0_ext = "",
  model_form_mu0_rct = "",
  model_form_mu1_rct = ""
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

\item{model_form_mu0_ext}{}

\item{model_form_mu0_rct}{}

\item{model_form_mu1_rct}{}
}
\value{
tau and standard deviation
}
\description{
Difference in difference + outcome regression + external control borrowing
}
\examples{
model_form_mu = c("y1 ~ x1 + x2", 
"y2 ~ x1 + x2",
"y3 ~ x1 + x2",
"y4 ~ x1 + x2")
res1 = DID_EC_OR(outcome = Y, trial_status = S, treatment = A, covariates = X,
                 long_term_marker = c(F, T, T, T), 
                 model_form_mu0_ext = model_form_mu, 
                 model_form_mu0_rct = model_form_mu, 
                 model_form_mu1_rct = model_form_mu)

}
