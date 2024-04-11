#' Analysis OLE class
#'
#' @slot method_obj method_OLE_obj
#' @slot T_cross 
#'
#' @include method_DID_class.R
#' @include method_SCM_class.R
#' @include method_class.R
#' @include analysis_class.R
#' @export setup_analysis_OLE
#'
#' @examples
#' 
#' bootstrap_obj = setup_bootstrap(
#'   replicates = 2e3, 
#'   bootstrap_CI_type = "perc" 
#' )
#' 
#' method_DID_obj = setup_method_DID(
#'   method_name = "IPW",
#'   bootstrap_flag = T,
#'   bootstrap_obj = bootstrap_obj,
#'   model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5",
#'   model_form_piA = "A ~ x1 + x2 + x3 + x4 + x5")
#' 
#' analysis_OLE_obj = setup_analysis_OLE( 
#'   data = SyntheticData, 
#'   trial_status_col_name = "S", 
#'   treatment_col_name = "A", 
#'   outcome_col_name = c("y1", "y2", "y3", "y4"), 
#'   covariates_col_name = c("x1", "x2", "x3", "x4", "x5"), 
#'   T_cross = 2, 
#'   method_OLE_obj = method_DID_obj)

.analysis_OLE_obj = setClass(
  "analysis_OLE_obj",
  contains = "analysis_obj",
  slots = c(
    method_obj = "method_OLE_obj",
    T_cross = "numeric"
  )
)

## TODO: modify the show method
setMethod(
  f = "show",
  signature = "analysis_OLE_obj",
  definition = function(object) {
    full_data = object@data
    print(full_data)
    cat("Using method: ", object@method_obj@method_name)
  }
)

#' Construct an analysis_OLE object
#'
#' @param method_obj 
#' @param data 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param outcome_col_name 
#' @param covariates_col_name 
#' @param T_cross 
#' @param alpha 
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
setup_analysis_OLE = function(data, trial_status_col_name, treatment_col_name, 
                          outcome_col_name, covariates_col_name, method_OLE_obj,
                          T_cross, alpha = 0.05){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  # length of long_term_marker the same as outcome dimension
  # if long_term_flag = F, then long_term_marker should all be F
  
  analysis_OLE_obj = .analysis_OLE_obj(
    data = data,
    covariates_col_name = covariates_col_name,
    outcome_col_name = outcome_col_name, 
    treatment_col_name = treatment_col_name, 
    trial_status_col_name = trial_status_col_name,
    method_obj = method_OLE_obj,
    T_cross = T_cross,
    alpha = alpha
  )
  
}