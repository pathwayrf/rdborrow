% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_DID_class.R
\docType{class}
\name{method_DID_obj-class}
\alias{method_DID_obj-class}
\alias{.method_DID_obj}
\title{method class}
\description{
method class
}
\section{Slots}{

\describe{
\item{\code{model_form_mu0_ext}}{}

\item{\code{model_form_mu0_rct}}{}

\item{\code{model_form_mu1_rct}}{}

\item{\code{model_form_piS}}{character.}

\item{\code{bootstrap_flag}}{}

\item{\code{model_form_piA}}{}

\item{\code{bootstrap_obj}}{}
}}

\examples{
method_DID_obj = setup_method_DID(
   method_name = "IPW",
   bootstrap_flag = T,
   bootstrap_obj = bootstrap_obj,
   model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5",
   model_form_piA = "A ~ x1 + x2 + x3 + x4 + x5")

}
