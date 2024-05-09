#' Difference in difference + AIPW + external control borrowing
#'
#' @param model_form_piS 
#' @param model_form_piA 
#' @param model_form_mu0_ext 
#' @param data 
#' @param indices 
#' @param outcome_col_name 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param covariates_col_name 
#' @param T_cross 
#'
#' @return tau and standard deviation
#' @export
#'
#' @examples 
#' DID_EC_AIPW_bootstrap = function(data,
#'   indices,
#'   outcome_col_name, 
#'   trial_status_col_name, 
#'   treatment_col_name, 
#'   covariates_col_name,
#'   T_cross,
#'   model_form_piS = "",
#'   model_form_piA = "",
#'   model_form_mu0_ext = ""
#' )
DID_EC_AIPW_bootstrap = function(data,
                       indices,
                       outcome_col_name, 
                       trial_status_col_name, 
                       treatment_col_name, 
                       covariates_col_name,
                       T_cross,
                       model_form_piS = "",
                       model_form_piA = "",
                       model_form_mu0_ext = ""
                     ){
  
  Y = subset(data[indices, ], select = outcome_col_name)
  S = subset(data[indices, ], select = trial_status_col_name)
  A = subset(data[indices, ], select = treatment_col_name)
  X = subset(data[indices, ], select = covariates_col_name)
  
  df = data[indices, ]
  N = nrow(df)
  T_follow = ncol(Y)
  T_pc = T_cross# placebo-control period
  n = sum(df$S) # RCT sample size
  m = sum(1 - df$S) # external control sample size
  pi.S = n/N
  
  # 
  piS = glm(as.formula(model_form_piS), data = df, family = "binomial")
  suppressWarnings({piSX = predict(piS, newdata = filter(df), type = "response")})
  if (model_form_piA == ""){
    piAX = sum(A)/n
  }else{
    piA = glm(as.formula(model_form_piA), data = filter(df, S == 1), family = "binomial")
    suppressWarnings({piAX = predict(piA, newdata = filter(df), type = "response")})
  }
  
  # predict Y0 from outcome regression models
  model_list_ext = lapply(1:T_follow, function(x){
    assign(paste0("m.ext", x), lm(as.formula(model_form_mu0_ext[x]), data = filter(df, S==0)))
  })
  suppressWarnings({Y0 = data.frame(sapply(1:T_follow, function(x){predict(model_list_ext[[x]], newdata = filter(df))}))})
  colnames(Y0) = paste0("y", 1:T_follow, "_0")
  # for residual
  Yr = Y - Y0 
  colnames(Yr) = paste0("y", 1:T_follow, "_r")
  
  temp = df %>% cbind(., Y0, Yr) %>%
    mutate(piAX = piAX, #sum(A)/sum(S)
           piSX = piSX,
           rx = piSX * (1 - pi.S)/(1 - piSX)/pi.S) %>%
    mutate(w11 = 1/piAX, w10 = 1/(1 - piAX), w00 = rx)
  
  # create outcomes
  Ys = as.matrix(Yr)
  potential = (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11) +
                 temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10) +
                 (1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00)) * Ys
  
  
  mu_S1A1 = colSums(data.frame(potential[temp$S == 1 & temp$A == 1, (T_pc+1):T_follow]))
  mu_S0A0 = colSums(data.frame(potential[temp$S == 0, (T_pc+1):T_follow]))
  bias = sum(rowMeans(data.frame(potential[temp$S==1 & temp$A==0, 1:T_pc]))) - 
    sum(rowMeans(data.frame(potential[temp$S==0, 1:T_pc])))
  
  
  tau = mu_S1A1 - mu_S0A0 - bias
  
  names(tau) = paste0("tau", (T_pc+1):T_follow)
  
  return(tau)
  
}

