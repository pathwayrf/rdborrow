#' Simulation class
#'
#' @slot covariates_col_name 
#' @slot outcome_col_name 
#' @slot treatment_col_name 
#' @slot trial_status_col_name 
#' @slot alpha 
#' @slot method_obj_list 
#' @slot method_description 
#'
#' @include method_class.R
#' @export setup_simulation
#'
#' 

 
.simulation_obj = setClass(
  "simulation_obj",
  slots = c(
    covariates_col_name = "character",
    outcome_col_name = "character",
    treatment_col_name = "character",
    trial_status_col_name = "character",
    method_obj_list = "list", # should be a list of method_obj object
    alpha = "numeric",
    method_description = "character"
  ),
  prototype = list(
    covariates_col_name = ""
  )
)



#' Simulation for primary analysis
#'
#' @slot data_matrix_list_null 
#' @slot data_matrix_list_alt 
#' @slot true_effect 
#' @slot alt_effect 
#'
#' @return a simulation object
#' @export
#'

.simulation_primary_obj = setClass(
  "simulation_primary_obj",
  contains = "simulation_obj",
  slots = c(
    data_matrix_list_null = "list",
    data_matrix_list_alt = "list",
    true_effect = "numeric",
    alt_effect = "numeric"
  ),
  prototype = list(
    data_matrix_list_alt = list(),
    alt_effect = numeric(0)
  )
)



#' Simulation for OLE study
#'
#' @slot T_cross numeric. 
#'
#' @return a simulation object for OLE phase
#' @export
#'
#' 
#' 
.simulation_OLE_obj = setClass(
  "simulation_OLE_obj",
  contains = "simulation_obj",
  slots = c(
    data_matrix_list = "list",
    true_effect = "numeric",
    T_cross = "numeric"
  )
)

## to be prepared for sanity check:
# data list: should be a list of data frames
# method_obj_list: should be a list of method_obj objects

## TODO: modify the show method
# setMethod(
#   f = "show",
#   signature = "simulation_obj",
#   definition = function(object) {
#     full_data = data
#     print(full_data)
#     cat("Using method: ", object@method)
#   }
# )

#' Construct an simulation object
#'
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param outcome_col_name 
#' @param covariates_col_name 
#' @param method_obj_list 
#' @param alpha 
#' @param method_description 
#'
#' @return An simulation object 
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
setup_simulation = function(trial_status_col_name, 
                            treatment_col_name, 
                            outcome_col_name, 
                            covariates_col_name, 
                            method_obj_list,
                            method_description,
                            alpha = 0.05){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  
  simulation_obj = .simulation_obj(
    covariates_col_name = covariates_col_name,
    outcome_col_name = outcome_col_name, 
    treatment_col_name = treatment_col_name, 
    trial_status_col_name = trial_status_col_name,
    method_obj_list = method_obj_list,
    alpha = alpha,
    method_description = method_description
  )
  
}


#' Construct a simulation object for primary analysis
#'
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param outcome_col_name 
#' @param covariates_col_name 
#' @param method_obj_list 
#' @param true_effect 
#' @param method_description 
#' @param alpha 
#' @param data_matrix_list_null 
#' @param data_matrix_list_alt 
#' @param alt_effect 
#'
#' @return return a simulation object for primary analysis
#' @export
#'
#' @examples
#' simulation_primary_obj = setup_simulation_primary(
#'  data_matrix_list_null = data_matrix_list_null,  # two scenarios
#'  data_matrix_list_alt = data_matrix_list_alt,
#'  trial_status_col_name = trial_status_col_name, 
#'  treatment_col_name = treatment_col_name, 
#'  outcome_col_name = outcome_col_name, 
#'  covariates_col_name = covariates_col_name, 
#'  method_obj_list = method_obj_list, 
#'  true_effect = true_effect,  
#'  alt_effect = alt_effect,
#'  alpha = alpha,
#'  method_description = c("IPW, optimal weight", 
#'                         "AIPW, optimal weight", 
#'                         "IPW, zero weight",
#'                         "AIPW, zero weight"))
#' 
setup_simulation_primary = function(data_matrix_list_null, 
                                    trial_status_col_name,
                                    treatment_col_name, 
                                    outcome_col_name, 
                                    covariates_col_name, 
                                    method_obj_list, 
                                    true_effect, 
                                    method_description,
                                    data_matrix_list_alt = list(),
                                    alt_effect = numeric(0),
                                    alpha = 0.05){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  
  simulation_primary_obj = .simulation_primary_obj(
    data_matrix_list_null = data_matrix_list_null,
    data_matrix_list_alt = data_matrix_list_alt,
    covariates_col_name = covariates_col_name,
    outcome_col_name = outcome_col_name, 
    treatment_col_name = treatment_col_name, 
    trial_status_col_name = trial_status_col_name,
    method_obj_list = method_obj_list,
    alpha = alpha,
    true_effect = true_effect,
    alt_effect = alt_effect,
    method_description = method_description
  )
  
}


#' Construct a simulation object for OLE analysis
#'
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param outcome_col_name 
#' @param covariates_col_name 
#' @param method_obj_list 
#' @param T_cross 
#' @param true_effect 
#' @param method_description 
#' @param alpha 
#' @param data_matrix_list
#'
#' @return a simulation object for OLE phase
#' @export
#'
#' @examples
#' 
#' simulation_OLE_obj = setup_simulation_OLE(
#'   data_matrix_list = data_matrix_list,  # two scenarios
#'   trial_status_col_name = trial_status_col_name, 
#'   treatment_col_name = treatment_col_name, 
#'   outcome_col_name = outcome_col_name, 
#'   covariates_col_name = covariates_col_name, 
#'   method_obj_list = method_obj_list, 
#'   true_effect = true_effect_long,  
#'   T_cross = 2,
#'   alpha = alpha,
#'   method_description = c("IPW, DID", 
#'                          "AIPW, DID", 
#'                          "OR, DID"))
#' 
#' 
#' 
setup_simulation_OLE = function(data_matrix_list, 
                                      trial_status_col_name, 
                                      treatment_col_name, 
                                      outcome_col_name, 
                                      covariates_col_name, 
                                      method_obj_list, 
                                      T_cross, 
                                      true_effect, 
                                      method_description,
                                      alpha = 0.05){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  
  simulation_OLE_obj = .simulation_OLE_obj(
    data_matrix_list = data_matrix_list, 
    covariates_col_name = covariates_col_name, 
    outcome_col_name = outcome_col_name, 
    treatment_col_name = treatment_col_name, 
    trial_status_col_name = trial_status_col_name, 
    method_obj_list = method_obj_list, 
    alpha = alpha, 
    T_cross = T_cross,
    true_effect = true_effect,
    method_description = method_description
  )
  
}

