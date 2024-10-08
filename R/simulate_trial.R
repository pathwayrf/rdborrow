#' Simulate trials
#'
#' @param X_int 
#' @param X_ext 
#' @param outcome_model_specs 
#' @param num_treated 
#' @param OLE_flag 
#' @param T_cross 
#'
#' @return a data frame for the simulated data
#' @export
#'
#' @examples
#' 
#' Data = simulate_trial(X_int, 
#'  X_ext, 
#'  num_treated = 150, 
#'  OLE_flag = T, 
#'  T_cross = 2,  
#'  outcome_model_specs) 
#' 
simulate_trial = function(X_int, X_ext, num_treated, OLE_flag, T_cross, outcome_model_specs){
  # TODO: sanity check for dimension of num_treated, T_cross, outcome_model_specs
  
  # ===== calculate population quantites =====
  n_int = nrow(X_int)
  n_ext = nrow(X_ext)
  T_follow = length(outcome_model_specs)
  
  # ===== Add treatment allocation indicator ===== 
  ## for internal group, A ~ Bernoulli(2/3) 
  # A_int = rbinom(n_int, 1, prob_treated)
  A_int = sample(c(rep(1, num_treated), rep(0, n_int - num_treated)))

  ## for external group, A = 0
  A_ext = rep(0, n_ext) 
  
  # ===== Add trial participation status indicator =====
  ## for internal group, S = 1
  S_int = rep(1, n_int) 
  
  ## for external group, S = 0 
  S_ext = rep(0, n_ext) 
  
  # ===== Simulate outcomes =====
  ## For internal data S = 1:
  Y_int = simulate_outcome_from_model(OLE_flag = OLE_flag,
                                      T_cross = T_cross, # number of followup time points
                                      X = X_int,    # covariates are specified by X_int
                                      A = A_int,    # treatment specified by A_int
                                      outcome_model_specs) # specify outcome model
  
  ## For external data S = 0:
  Y_ext = simulate_outcome_from_model(OLE_flag = F,
                                      T_cross = T_cross, # number of followup time points
                                      X = X_ext,    # covariates are specified by X_ext
                                      A = A_ext,    # treatment is specified by A_ext
                                      outcome_model_specs) # specify outcome model
  
  # ===== Combine two different studies =====
  if(OLE_flag){
    ## combine pieces for the internal group S = 1
    Data_int = data.frame(X_int, A = A_int, S = S_int, T_cross = rep(T_cross, n_int), Y_int)
    
    ## combine pieces for the external group S = 0
    Data_ext = data.frame(X_ext, A = A_ext, S = S_ext, T_cross = rep(T_cross, n_ext), Y_ext)
  }else{
    ## combine pieces for the internal group S = 1
    Data_int = data.frame(X_int, A = A_int, S = S_int, Y_int)
    
    ## combine pieces for the external group S = 0
    Data_ext = data.frame(X_ext, A = A_ext, S = S_ext, Y_ext)
  }
  
  
  ## combine internal and external data
  Data = rbind(Data_int, Data_ext)
  
  rownames(Data) = as.character(1:(n_int + n_ext))
  Data
  
}