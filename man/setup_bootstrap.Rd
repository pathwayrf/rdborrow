% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_class.R
\name{setup_bootstrap}
\alias{setup_bootstrap}
\title{Construct a bootstrap object}
\usage{
setup_bootstrap(replicates = 500, bootstrap_CI_type = "bca")
}
\arguments{
\item{replicates}{}

\item{bootstrap_CI_type}{}
}
\value{
An bootstrap object
}
\description{
Construct a bootstrap object
}
\examples{
analysis_obj = setup_analysis(trial_status_col_name = S, 
   treatment_col_name = A, 
   outcome_col_name = Y, 
   covariates_col_name = X, 
   method = method_obj)


}
