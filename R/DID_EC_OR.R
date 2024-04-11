#' Difference in difference + outcome regression + external control borrowing
#'
#' @param model_form_mu0_ext 
#' @param model_form_mu0_rct 
#' @param model_form_mu1_rct 
#' @param data 
#' @param outcome_col_name 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param covariates_col_name 
#' @param Bootstrap 
#' @param R 
#' @param bootstrap_CI_type 
#' @param alpha 
#' @param T_cross 
#'
#' @include DID_EC_OR_bootstrap.R
#' @return tau and standard deviation
#' @export
#'
#' @examples 
#' model_form_mu = c("y1 ~ x1 + x2", 
#' "y2 ~ x1 + x2",
#' "y3 ~ x1 + x2",
#' "y4 ~ x1 + x2")
#' res1 = DID_EC_OR(outcome = Y, trial_status = S, treatment = A, covariates = X,
#'                  long_term_marker = c(F, T, T, T), 
#'                  model_form_mu0_ext = model_form_mu, 
#'                  model_form_mu0_rct = model_form_mu, 
#'                  model_form_mu1_rct = model_form_mu)
#' 
DID_EC_OR = function(data,
                     outcome_col_name, 
                     trial_status_col_name, 
                     treatment_col_name, 
                     covariates_col_name,
                     T_cross,
                     model_form_mu0_ext = "",
                     model_form_mu0_rct = "",
                     model_form_mu1_rct = "",
                     Bootstrap = F,
                     R = 5e2,
                     bootstrap_CI_type = "bca", 
                     alpha = 0.05,
                     quiet = TRUE
                     ){
  
  Y = subset(data, select = outcome_col_name)
  S = subset(data, select = trial_status_col_name)
  A = subset(data, select = treatment_col_name)
  X = subset(data, select = covariates_col_name)
  
  df = data
  N = nrow(df)
  T_follow = ncol(Y)
  T_pc = T_cross # placebo-control period
  n = sum(df$S) # RCT sample size
  m = sum(1 - df$S) # external control sample size
  pi.S = n/N
  
  
  
  # external outcome model
  model_list_ext = lapply(1:T_follow, function(x){
    assign(paste0("m.ext", x), lm(as.formula(model_form_mu0_ext[x]), data = filter(df, S==0)))
  })
  
  # list2env(setNames(model.list, c("m.ext1","m.ext2","m.ext3","m.ext4")), envir = .GlobalEnv)
  
  # rct outcome model
  model_list_rct_pc = lapply(1:T_pc, function(x){
    assign(paste0("m.rct", x), lm(as.formula(model_form_mu0_rct[x]), data = filter(df, S==1 & A==0)))
  })
  
  
  model_list_rct_cr = lapply((T_pc+1):T_follow, function(x){
    assign(paste0("m.rct", x), lm(as.formula(model_form_mu1_rct[x]), data = filter(df, S==1 & A==1)))
  })
  
  model_list_rct = c(model_list_rct_pc, model_list_rct_cr)
  
  
  # predicted value for external subjects
  mu_S0A0 = do.call(cbind,
                    lapply(
                      1:T_follow,
                      function(x){
                        predict(model_list_ext[[x]], newdata = filter(df, S==1))
                      }
                    ))
  colnames(mu_S0A0) = paste0("mu_S0A0_", 1:T_follow)
  
  # S=1
  # first T_pc cols for mu_S1A0 for placebo-controlled period,
  # rest cols for mu_S1A1 for crossover stage
  mu_S1 = do.call(cbind,
                    lapply(
                      1:T_follow,
                      function(x){
                        predict(model_list_rct[[x]], newdata = filter(df, S==1))
                      }
                    ))
  colnames(mu_S1) = paste0("mu_S1_", 1:T_follow)
  
  avg_S0A0 = colMeans(mu_S0A0)
  avg_S1 = colMeans(mu_S1)
  
  tau = (avg_S1[(T_pc+1):T_follow] - mean(avg_S1[1:T_pc])) - 
    (avg_S0A0[(T_pc+1):T_follow] - mean(avg_S0A0[1:T_pc]))
  
  names(tau) = paste0("tau", (T_pc+1):T_follow)

  # summarize results
  
  cutoff = qnorm(1-alpha/2, lower.tail = T)
  
  if(Bootstrap){
    Group_ID = df %>% group_by(S, A) %>% mutate(group_id = cur_group_id())
    Group_ID = Group_ID$group_id
    
    boot.ci.type = switch (bootstrap_CI_type,
                           norm = "normal",
                           bca = "bca",
                           stud = "student",
                           perc = "percent",
                           basic = "basic"
    )
    
    boot.out <- boot(data = df, 
                     statistic = DID_EC_OR_bootstrap, 
                     outcome_col_name = outcome_col_name, 
                     trial_status_col_name = trial_status_col_name, 
                     treatment_col_name = treatment_col_name, 
                     covariates_col_name = covariates_col_name,
                     T_cross = T_cross,
                     model_form_mu0_ext = model_form_mu0_ext,
                     model_form_mu0_rct = model_form_mu0_rct,
                     model_form_mu1_rct = model_form_mu1_rct,
                     R = R, 
                     strata = Group_ID)
    
    lower_CI_boot = sapply(1:(T_follow - T_pc), 
                           function(x) {
                             ci = boot.ci(boot.out, conf = 1 - alpha, type = bootstrap_CI_type, index = x)
                             ci[[boot.ci.type]][4]
                           })
    
    upper_CI_boot = sapply(1:(T_follow - T_pc), 
                           function(x) {
                             ci = boot.ci(boot.out, conf = 1 - alpha, type = bootstrap_CI_type, index = x)
                             ci[[boot.ci.type]][5]
                           })
    
    sd.boot = sqrt(diag(var(boot.out$t)))
    
    results = data.frame(
      point_estimates = tau,
      lower_CI_boot = lower_CI_boot,
      upper_CI_boot = upper_CI_boot
    )
    
  }else{
    stop("No other inference methods defined!")
  }
  
  return(results)
  
}

