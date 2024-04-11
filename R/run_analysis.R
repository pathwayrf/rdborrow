#' run analysis
#' 
#' @param analysis_obj 
#' @param Bootstrap 
#' 
#' @return  a list containing: tau (effect size), sd.tau (standard deviation), wt (weight)
#' @include EC_IPW_OPT.R
#' @include EC_AIPW_OPT.R
#' @include DID_EC_IPW.R
#' @include DID_EC_AIPW.R
#' @include DID_EC_OR.R
#' @include SCM.R
#' 
#' @export
#' 
#' @examples 
#' run_analysis(analysis_obj, Bootstrap = F)

run_analysis = function(analysis_obj, quiet = TRUE){
  
  # sanity check
  ## TODO: 
  # correct initialization of objects
  # correct dimension compatible
  # 
  
  data = analysis_obj@data
  outcome_col_name = analysis_obj@outcome_col_name
  trial_status_col_name = analysis_obj@trial_status_col_name
  treatment_col_name = analysis_obj@treatment_col_name
  covariates_col_name = analysis_obj@covariates_col_name
  alpha = analysis_obj@alpha
  method = analysis_obj@method_obj
  Bootstrap = method@bootstrap_obj
  
  bootstrap_flag = method@bootstrap_flag
  R = Bootstrap@replicates
  bootstrap_CI_type = Bootstrap@bootstrap_CI_type
  
  method_type = is(method)[1]
  if (is(analysis_obj)[1] == "analysis_OLE_obj"){
    T_cross = analysis_obj@T_cross
  }
  
  if (method_type == "method_weighting_obj"){
    if (!quiet){
      cat("Estimating causal effects for primary analysis... \n")
    }
    
    wt = method@wt
    optimal_weight_flag = method@optimal_weight_flag
    
    name = method@method_name
    model_form_piS = method@model_form_piS
    model_form_piA = method@model_form_piA
    model_form_mu0_ext = method@model_form_mu0_ext
    model_form_mu0_rct = method@model_form_mu0_rct
    model_form_mu1_rct = method@model_form_mu1_rct
    
    
    if (name == "IPW"){
      res = EC_IPW_OPT(data = data,
                       outcome_col_name = outcome_col_name, 
                       trial_status_col_name = trial_status_col_name, 
                       treatment_col_name = treatment_col_name, 
                       covariates_col_name = covariates_col_name, 
                       model_form_piS = model_form_piS,
                       wt = wt, 
                       optimal_weight_flag = optimal_weight_flag,
                       bootstrap_flag, R, bootstrap_CI_type, alpha = alpha, quiet = quiet)
    }else if (name == "AIPW"){
      res = EC_AIPW_OPT(data = data,
                        outcome_col_name = outcome_col_name, 
                        trial_status_col_name = trial_status_col_name, 
                        treatment_col_name = treatment_col_name, 
                        covariates_col_name = covariates_col_name, 
                        model_form_piS = model_form_piS,
                        model_form_mu0_ext = model_form_mu0_ext,
                        wt = wt, 
                        optimal_weight_flag = optimal_weight_flag,
                        bootstrap_flag, R, bootstrap_CI_type, alpha = alpha, quiet = quiet)
    }else {
      stop("No such method is defined!")
    }
  }else if (method_type == "method_DID_obj") {
    ## TODO: implement OLE analysis
    if (!quiet){
      cat("Estimating long term causal effects using DID methods... \n")
    }
    
    name = method@method_name
    model_form_piS = method@model_form_piS
    model_form_piA = method@model_form_piA
    model_form_mu0_ext = method@model_form_mu0_ext
    model_form_mu0_rct = method@model_form_mu0_rct
    model_form_mu1_rct = method@model_form_mu1_rct
    
    
    if (name == "IPW"){
      res = DID_EC_IPW(data = data,
                       outcome_col_name = outcome_col_name, 
                       trial_status_col_name = trial_status_col_name, 
                       treatment_col_name = treatment_col_name, 
                       covariates_col_name = covariates_col_name,
                       T_cross = T_cross,
                       model_form_piS = model_form_piS,
                       model_form_piA = model_form_piA,
                       Bootstrap = bootstrap_flag,
                       R = R,
                       bootstrap_CI_type = bootstrap_CI_type, 
                       alpha = alpha, quiet = quiet)
    }else if (name == "AIPW"){
      res = DID_EC_AIPW(data = data,
                        outcome_col_name = outcome_col_name, 
                        trial_status_col_name = trial_status_col_name, 
                        treatment_col_name = treatment_col_name, 
                        covariates_col_name = covariates_col_name,
                        T_cross = T_cross,
                        model_form_piS = model_form_piS,
                        model_form_piA = model_form_piA,
                        model_form_mu0_ext = model_form_mu0_ext,
                        Bootstrap = bootstrap_flag,
                        R = R,
                        bootstrap_CI_type = bootstrap_CI_type, 
                        alpha = alpha, quiet = quiet)
    }else if (name == "OR"){
      res = DID_EC_OR(data = data,
                      outcome_col_name = outcome_col_name, 
                      trial_status_col_name = trial_status_col_name, 
                      treatment_col_name = treatment_col_name, 
                      covariates_col_name = covariates_col_name,
                      T_cross = T_cross,
                      model_form_mu0_ext = model_form_mu0_ext,
                      model_form_mu0_rct = model_form_mu0_rct,
                      model_form_mu1_rct = model_form_mu1_rct,
                      Bootstrap = bootstrap_flag,
                      R = R,
                      bootstrap_CI_type = bootstrap_CI_type, 
                      alpha = alpha, quiet = quiet)
    }
    else {
      stop("No such method is defined!")
    }
  }else if (method_type == "method_SCM_obj"){
    if(!quiet){
      cat("Estimating long term causal effects using synthetic control methods... \n")
    }
    
    name = method@method_name
    lambda.min = method@lambda.min
    lambda.max = method@lambda.max
    nlambda = method@nlambda
    parallel = method@parallel
    ncpus = method@ncpus
    
    if (name == "SCM"){
      res = SCM(data = data, 
                trial_status_col_name = trial_status_col_name, 
                treatment_col_name = treatment_col_name, 
                outcome_col_name = outcome_col_name,  
                covariates_col_name = covariates_col_name,  
                T_cross = T_cross, 
                Bootstrap = Bootstrap, 
                R = R, 
                bootstrap_CI_type = bootstrap_CI_type, 
                alpha = alpha, 
                lambda.min = lambda.min, 
                lambda.max = lambda.max, 
                nlambda = nlambda, 
                parallel = parallel, 
                ncpus = ncpus, quiet = quiet)
    } else{
      stop("No such method is defined!")
    }
  }else {
    stop("No such method type is defined!")
  }
  
  return(res)
  
}