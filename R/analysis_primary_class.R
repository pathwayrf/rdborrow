#' Analysis class
#'
#' @slot method_obj 
#'
#' @include analysis_class.R
#' @include method_weighting_class.R
#' @export setup_analysis_primary
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
#' analysis_primary_obj = setup_analysis_primary(
#'   data = SyntheticData,
#'   trial_status_col_name = "S", 
#'   treatment_col_name = "A", 
#'   outcome_col_name = c("y1", "y2"), 
#'   covariates_col_name = c("x1", "x2", "x3", "x4", "x5"), 
#'   method_weighting_obj = method_weighting_obj)
#' 

.analysis_primary_obj = setClass(
  "analysis_primary_obj",
  contains = "analysis_obj",
  slots = c(
    method_obj = "method_primary_obj"
  )
)

## TODO: modify the show method
setMethod(
  f = "show",
  signature = "analysis_primary_obj",
  definition = function(object) {
    full_data = data
    print(full_data)
    cat("Using method: ", object@method)
  }
)

#' Construct an analysis_OLE object
#'
#' @param data 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param outcome_col_name 
#' @param covariates_col_name 
#' @param alpha 
#' @param method_weighting_obj 
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
setup_analysis_primary = function(data, trial_status_col_name, treatment_col_name, 
                          outcome_col_name, covariates_col_name, method_weighting_obj,
                          alpha = 0.05){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  # length of long_term_marker the same as outcome dimension
  # if long_term_flag = F, then long_term_marker should all be F
  
  analysis_primary_obj = .analysis_primary_obj(
    data = data,
    covariates_col_name = covariates_col_name,
    outcome_col_name = outcome_col_name, 
    treatment_col_name = treatment_col_name, 
    trial_status_col_name = trial_status_col_name,
    method_obj = method_weighting_obj,
    alpha = alpha
  )
  
}