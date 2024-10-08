% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCM.R
\name{SCM}
\alias{SCM}
\title{Implement the Synthetic Control Method}
\usage{
SCM(
  data,
  outcome_col_name,
  trial_status_col_name,
  treatment_col_name,
  covariates_col_name,
  T_cross,
  Bootstrap = T,
  R = 100,
  bootstrap_CI_type = "bca",
  alpha = 0.05,
  lambda.min = 0,
  lambda.max = 0.1,
  nlambda = 2,
  parallel = "no",
  ncpus = 1,
  quiet = TRUE
)
}
\arguments{
\item{data}{A data frame}

\item{outcome_col_name}{}

\item{trial_status_col_name}{}

\item{treatment_col_name}{}

\item{covariates_col_name}{}

\item{T_cross}{}

\item{Bootstrap}{}

\item{R}{}

\item{bootstrap_CI_type}{}

\item{alpha}{}

\item{lambda.min}{}

\item{lambda.max}{}

\item{nlambda}{}

\item{parallel}{}

\item{ncpus}{}
}
\value{
A list contains: estimated ATE, SE, weight used, SE by Bootstrap and a 95% confidence interval for primary endpoint (only when Bootstrap=TRUE)
}
\description{
SCM() is the main function calculates the estimated ATE by SC method and Bootstrap CI, it calls subject_SC() and lambdacv().
}
\examples{
This is part of the run_analysis() function

 
}
