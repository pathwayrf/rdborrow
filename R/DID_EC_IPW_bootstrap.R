#' Difference in difference + IPW + external control borrowing
#'
#' @param model_form_piS 
#' @param data 
#' @param indices 
#' @param outcome_col_name 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param covariates_col_name 
#' @param T_cross 
#' @param model_form_piA 
#'
#' @return tau and standard deviation
#' @export
#'
#' @examples This is part of the run_analysis() function
#' 
DID_EC_IPW_bootstrap = function(data,
                      indices,
                      outcome_col_name, 
                      trial_status_col_name, 
                      treatment_col_name, 
                      covariates_col_name,
                      T_cross,
                      model_form_piS = "",
                      model_form_piA = ""
                     ){
  
  Y = subset(data[indices, ], select = outcome_col_name)
  S = subset(data[indices, ], select = trial_status_col_name)
  A = subset(data[indices, ], select = treatment_col_name)
  X = subset(data[indices, ], select = covariates_col_name)
  
  df = data[indices, ]
  N = nrow(df)
  T_follow = ncol(Y)
  T_pc = T_cross # placebo-control period
  n = sum(df$S) # RCT sample size
  m = sum(1 - df$S) # external control sample size
  pi.S = n/N
  
  # 
  piS = glm(as.formula(model_form_piS), data = df, family = "binomial")
  piSX = predict(piS, newdata = filter(df), type = "response")
  if (model_form_piA == ""){
    piAX = sum(A)/n
  }else{
    piA = glm(as.formula(model_form_piA), data = filter(df, S == 1), family = "binomial")
    piAX = predict(piA, newdata = filter(df), type = "response")
  }
  
  
  temp = df %>% 
    mutate(piAX = piAX, #sum(A)/sum(S)
           piSX = piSX,
           rx = piSX * (1 - pi.S)/(1 - piSX)/pi.S) %>%
    mutate(w11 = 1/piAX, w10 = 1/(1 - piAX), w00 = rx)

  
  # create outcomes
  Ys = as.matrix(Y)
  potential = (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11) +
                 temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10) +
                 (1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00)) * Ys
  
  
  mu_S1A1 = colSums(data.frame(potential[temp$S == 1 & temp$A == 1, (T_pc+1):T_follow, drop = F]))
  mu_S0A0 = colSums(data.frame(potential[temp$S == 0, (T_pc+1):T_follow, drop = F]))
  bias = sum(rowMeans(data.frame(potential[temp$S==1 & temp$A==0, 1:T_pc, drop = F]))) - 
    sum(rowMeans(data.frame(potential[temp$S==0, 1:T_pc, drop = F])))
  
  
  tau = mu_S1A1 - mu_S0A0 - bias
  
  names(tau) = paste0("tau", (T_pc+1):T_follow)
  
  ####### Use Bootstrap for standard error and confidence intervals
  
  # boot.out <- boot(data=df, DiDboot,
  #                  setting = setting, 
  #                  method = method,
  #                  form_x = form_x,
  #                  form_x_omit = form_x_omit,
  #                  # parallel = "multicore",
  #                  # ncpus=4,
  #                  strata = df$S,
  #                  R = 1000)
  # 
  # 
  # return(list(c(tau[1], sd(boot.out$t[,1]), boot.ci(boot.out,index=1,type=c("perc"))$percent[4:5]),
  #             c(tau[2], sd(boot.out$t[,2]), boot.ci(boot.out,index=2,type=c("perc"))$percent[4:5])))
  
  return(tau)
  
}

