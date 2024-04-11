#' Run simulation from a simualtion_obj
#'
#' @param simulation_obj 
#' @param quiet 
#'
#' @return a simulation_report object
#' @export
#'
#' @examples
#' simulation_report = run_simulation(simulation_OLE_obj, quiet = F)
run_simulation = function(simulation_obj, quiet = TRUE){
  
  covariates_col_name = simulation_obj@covariates_col_name
  outcome_col_name = simulation_obj@outcome_col_name
  treatment_col_name = simulation_obj@treatment_col_name
  trial_status_col_name = simulation_obj@trial_status_col_name
  method_obj_list = simulation_obj@method_obj_list # should be a list of method_obj object
  alpha = simulation_obj@alpha
  
  nmethod = length(method_obj_list)
  power = numeric(0)
  
  record = list()
  
  # For null analysis
  
  if (is(simulation_obj)[1] == "simulation_primary_obj"){
    # for null analysis
    data_matrix_list_null = simulation_obj@data_matrix_list_null
    data_matrix_list_alt = simulation_obj@data_matrix_list_alt
    
    true_effect = simulation_obj@true_effect
    alt_effect = simulation_obj@alt_effect
    for (md_iter in 1:nmethod){
      method = method_obj_list[[md_iter]]
      record_df = data.frame()
      niter = length(data_matrix_list_null)
      for (iter in 1:niter){
        if (!quiet){
          cat("Null: ", "| method setting: ", md_iter, "| data setting: ", iter, "\n")
        }
        # create a primary analysis object
        analysis_primary_obj = setup_analysis_primary(
          data = data_matrix_list_null[[iter]],
          trial_status_col_name = trial_status_col_name, 
          treatment_col_name = treatment_col_name, 
          outcome_col_name = outcome_col_name, 
          covariates_col_name = covariates_col_name, 
          method_weighting_obj = method,
          alpha = alpha)
        
        res = run_analysis(analysis_primary_obj)
        rownames(res$results)[nrow(res$results)] = ""
        
        # only focus on the last time point
        # print(nrow(res$results))
        record_df = rbind(record_df, res$results[nrow(res$results), ])
      }
      record[[md_iter]] = record_df
    }
    
    method_description = simulation_obj@method_description
    bias = sapply(1:length(method_obj_list), function(x) { mean(record[[x]]$point_estimates - true_effect) })
    variance = sapply(1:length(method_obj_list), function(x) { var(record[[x]]$point_estimates) })
    mse = variance + bias^2
    coverage = sapply(1:length(method_obj_list), 
                      function(x){
                        if (("lower_CI_boot" %in%  colnames(record[[x]])) & ("upper_CI_boot" %in%  colnames(record[[x]]))){
                          return(sum((true_effect > record[[x]]$lower_CI_boot) & (true_effect < record[[x]]$upper_CI_boot))/length(data_matrix_list_null))
                        }else{
                          return(sum((true_effect > record[[x]]$lower_CI_normal) & (true_effect < record[[x]]$upper_CI_normal))/length(data_matrix_list_null))
                        }
                      })
    type_I_error = 1 - coverage
    
    # for power analysis
    if(length(data_matrix_list_alt) != 0){
      record_alt = list()
      niter = length(data_matrix_list_alt)
      for (md_iter in 1:nmethod){
        method = method_obj_list[[md_iter]]
        record_df_alt = data.frame()
        for (iter in 1:niter){
          if (!quiet){
            cat("Alternative: ", "| method setting: ", md_iter, "| data setting: ", iter, "\n")
          }
          # create a primary analysis object
          analysis_primary_obj = setup_analysis_primary(
            data = data_matrix_list_alt[[iter]],
            trial_status_col_name = trial_status_col_name, 
            treatment_col_name = treatment_col_name, 
            outcome_col_name = outcome_col_name, 
            covariates_col_name = covariates_col_name, 
            method_weighting_obj = method,
            alpha = alpha)
          
          res = run_analysis(analysis_primary_obj)
          rownames(res$results)[nrow(res$results)] = ""
          
          # only focus on the last time point
          # print(nrow(res$results))
          record_df_alt = rbind(record_df_alt, res$results[nrow(res$results), ])
        }
        record_alt[[md_iter]] = record_df_alt
      }
      

      names(record_alt) = paste0("method_candidate_", 1:length(method_obj_list))
      method_description = simulation_obj@method_description
      # bias = sapply(1:length(method_obj_list), function(x) { mean(record[[x]]$point_estimates - true_effect) })
      # variance = sapply(1:length(method_obj_list), function(x) { var(record[[x]]$point_estimates) })
      # mse = variance + bias^2
      power = 1 - sapply(1:length(method_obj_list), 
                         function(x){
                           if (("lower_CI_boot" %in%  colnames(record_alt[[x]])) & ("upper_CI_boot" %in%  colnames(record_alt[[x]]))){
                             return(sum((true_effect > record_alt[[x]]$lower_CI_boot) & (true_effect < record_alt[[x]]$upper_CI_boot))/length(data_matrix_list_alt))
                           }else{
                             return(sum((true_effect > record_alt[[x]]$lower_CI_normal) & (true_effect < record_alt[[x]]$upper_CI_normal))/length(data_matrix_list_alt))
                           }
                         })
    }
    
  }else if (is(simulation_obj)[1] == "simulation_OLE_obj"){
    data_matrix_list = simulation_obj@data_matrix_list
    niter = length(data_matrix_list)
    true_effect = simulation_obj@true_effect
    
    for (md_iter in 1:nmethod){
      method = method_obj_list[[md_iter]]
      record_df = data.frame()
      for (iter in 1:niter){
        if (!quiet){
          cat("Null: ", "| method setting: ", md_iter, "| data setting: ", iter, "\n")
        }
        T_cross = simulation_obj@T_cross
        # create a primary analysis object
        analysis_OLE_obj = setup_analysis_OLE(
          data = data_matrix_list[[iter]],
          trial_status_col_name = trial_status_col_name, 
          treatment_col_name = treatment_col_name, 
          outcome_col_name = outcome_col_name, 
          covariates_col_name = covariates_col_name, 
          method_OLE_obj = method,
          alpha = alpha,
          T_cross = T_cross)
        
        res = run_analysis(analysis_OLE_obj)
        # print(res)
        rownames(res)[nrow(res)] = ""
        
        # only focus on the last time point
        record_df = rbind(record_df, res[nrow(res), ])
      }
      record[[md_iter]] = record_df
    }
    
    method_description = simulation_obj@method_description
    bias = sapply(1:length(method_obj_list), function(x) { mean(record[[x]]$point_estimates - true_effect) })
    variance = sapply(1:length(method_obj_list), function(x) { var(record[[x]]$point_estimates) })
    mse = variance + bias^2
    coverage = sapply(1:length(method_obj_list), 
                      function(x){
                        if (("lower_CI_boot" %in%  colnames(record[[x]])) & ("upper_CI_boot" %in%  colnames(record[[x]]))){
                          return(sum((true_effect > record[[x]]$lower_CI_boot) & (true_effect < record[[x]]$upper_CI_boot))/length(data_matrix_list))
                        }else{
                          return(sum((true_effect > record[[x]]$lower_CI_normal) & (true_effect < record[[x]]$upper_CI_normal))/length(data_matrix_list))
                        }
                      })
    type_I_error = 1 - coverage
    
  }
  else{
    stop("No such objects are defined!")
  }
  
  # summary
  names(record) = paste0("method_candidate_", 1:length(method_obj_list))
  
  
  
  
  setup_simulation_report(method_description = method_description,
                          bias = bias,
                          variance = variance,
                          mse = mse,
                          coverage = coverage,
                          type_I_error = type_I_error,
                          power = power)
}