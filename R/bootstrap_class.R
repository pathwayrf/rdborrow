#' Bootstrap class
#'
#' @slot replicates 
#' @slot bootstrap_CI_type 
#'
#' @include method_class.R
#' @export setup_analysis
#'
#' @examples
#' bootstrap_obj = setup_bootstrap(
#'   replicates = 2e3,
#'   bootstrap_CI_type = "perc"
#' )
#' 
#' method_weighting_obj = setup_method_weighting(
#'   method_name = "AIPW",
#'   optimal_weight_flag = T, 
#'   wt = 0,
#'   bootstrap_flag = T,
#'   bootstrap_obj = bootstrap_obj,
#'   model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5",
#'   model_form_mu0_ext = c("y1 ~ x1 + x2 + x3 + x4 + x5",
#'                          "y2 ~ x1 + x2 + x3 + x4 + x5"))
#' 
#' analysis_primary_obj = setup_analysis_primary(
#'   data = SyntheticData,
#'   trial_status = "S", 
#'   treatment = "A", 
#'   outcome = c("y1", "y2"), 
#'   covariates = c("x1", "x2", "x3", "x4", "x5"), 
#'   method_weighting_obj = method_weighting_obj)
#'
#'
#' res = run_analysis(analysis_primary_obj)

.bootstrap_obj <- setClass(
  "bootstrap_obj",
  slots = c(
    replicates = "numeric",
    bootstrap_CI_type = "character"
  ),
    prototype = list(
    replicates = 5e2,
    bootstrap_CI_type = "bca"
  )
)

## TODO: modify the show method
setMethod(
  f = "show",
  signature = "bootstrap_obj",
  definition = function(object) {
#    full_data = data
#    print(full_data)
    # "norm","basic", "stud", "perc", "bca"
    boot.ci.type = switch (object@bootstrap_CI_type,
      norm = "normal approximation",
      bca = "bias-corrected",
      stud = "studentized",
      perc = "percentile",
      basic = "basic"
    )
    cat("Running Bootstrap: ", ifelse(object@bootstrap_flag, "Yes", "No"), "\n")
    cat("Number of Replicates: ", as.character(object@replicates), "\n")
    cat("Type of bootstrap confidence interval: ", boot.ci.type)
  }
)

#' Construct a bootstrap object
#'
#' @param replicates 
#' @param bootstrap_CI_type 
#'
#' @return An bootstrap object 
#' @export
#'
#' @examples
#' analysis_obj = setup_analysis(trial_status_col_name = S, 
#'    treatment_col_name = A, 
#'    outcome_col_name = Y, 
#'    covariates_col_name = X, 
#'    method = method_obj)
#' 
#' 
setup_bootstrap = function(replicates = 5e2,
                           bootstrap_CI_type = "bca"){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  # length of long_term_marker the same as outcome dimension
  # if long_term_flag = F, then long_term_marker should all be F
  
  bootstrap_obj = .bootstrap_obj(
    replicates = replicates,
    bootstrap_CI_type = bootstrap_CI_type
  )
  
}