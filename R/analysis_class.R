#' Analysis class
#'
#' @slot method_obj Method. 
#' @slot data 
#' @slot covariates_col_name 
#' @slot outcome_col_name 
#' @slot treatment_col_name 
#' @slot trial_status_col_name 
#' @slot alpha 
#'
#' @include method_class.R
#' @export setup_analysis
#'
#' @examples
#' 
#' method_weighting_obj = setup_method_weighting(
#'   method_name = "IPW",
#'   optimal_weight_flag = F, 
#'   wt = 0,
#'   model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5")
#' 
#' 
#' analysis_obj = setup_analysis(
#'   data = SyntheticData,
#'   trial_status_col_name = "S", 
#'   treatment_col_name = "A", 
#'   outcome_col_name = c("y1", "y2"), 
#'   covariates_col_name = c("x1", "x2", "x3", "x4", "x5"), 
#'   method_obj = method_weighting_obj)
#' 
 
.analysis_obj = setClass(
  "analysis_obj",
  slots = c(
    data = "data.frame",
    covariates_col_name = "character",
    outcome_col_name = "character",
    treatment_col_name = "character",
    trial_status_col_name = "character",
    method_obj = "method_obj", # change to method_obj
    alpha = "numeric"
  ),
  prototype = list(
    covariates_col_name = ""
  )
)

## TODO: modify the show method
setMethod(
  f = "show",
  signature = "analysis_obj",
  definition = function(object) {
    full_data = data
    print(full_data)
    cat("Using method: ", object@method)
  }
)

#' Construct an analysis object
#'
#' @param method_obj 
#' @param data 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param outcome_col_name 
#' @param covariates_col_name 
#'
#' @return An analysis object [`Analysis`][Analysis-class]
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
setup_analysis = function(data, trial_status_col_name, treatment_col_name, 
                          outcome_col_name, covariates_col_name, method_obj,
                          alpha = 0.05){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  
  analysis_obj = .analysis_obj(
    data = data,
    covariates_col_name = covariates_col_name,
    outcome_col_name = outcome_col_name, 
    treatment_col_name = treatment_col_name, 
    trial_status_col_name = trial_status_col_name,
    method_obj = method_obj,
    alpha = alpha
  )
  
}