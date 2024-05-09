#' Using AIPW with external borrowing
#'
#' @param model_form_piS 
#' @param model_form_mu0_ext 
#' @param optimal_weight_flag 
#' @param wt 
#' @param data 
#' @param outcome_col_name 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param covariates_col_name 
#'
#' @return a list containing: tau (effect size), sd.tau (standard deviation), wt (weight)
#' @export
#' @examples 
#' res = EC_AIPW_OPT(data = data,
#'   outcome_col_name = outcome_col_name, 
#'   trial_status_col_name = trial_status_col_name, 
#'   treatment_col_name = treatment_col_name, 
#'   covariates_col_name = covariates_col_name, 
#'   model_form_piS = model_form_piS,
#'   model_form_mu0 = model_form_mu0,
#'   wt = wt, 
#'   optimal_weight_flag = optimal_weight_flag,
#'   Bootstrap)
#' 

EC_AIPW_OPT_bootstrap = function(data,
                                 indices,
                                 outcome_col_name, 
                                 trial_status_col_name, 
                                 treatment_col_name, 
                                 covariates_col_name,
                                 model_form_piS = "",
                                 model_form_mu0_ext = "",
                                 optimal_weight_flag = F,
                                 wt = 0){
  Y = subset(data[indices, ], select = outcome_col_name)
  S = subset(data[indices, ], select = trial_status_col_name)
  A = subset(data[indices, ], select = treatment_col_name)
  X = subset(data[indices, ], select = covariates_col_name)
  
  df = data[indices, ]
  N = nrow(df)
  T_follow = ncol(Y)
  n = sum(df$S) # RCT sample size
  m = sum(1 - df$S) # external control sample size
  pi.S = n/N
  
  ## TODO: write running messages
  ## cat("running AIPW... \n")
  
  if((!optimal_weight_flag) && wt == 0){
    # estimate ATE
    ## TODO: why do use the true propensity score? 
    temp = df %>%
      filter(S==1) %>%
      mutate(`piA`=sum(A)/n) %>%
      mutate(w11 = `piA`, w10 = 1 - `piA`)
    
    ### create outcomes: obs by T
    Ys = as.matrix(Y[S==1, ])
    
    potential = (temp$A/temp$w11+(1-temp$A)/temp$w10)*Ys
    mu1 = colSums(potential[temp$A==1, ])/n
    mu0 = colSums(potential[temp$A==0, ])/n
    
    tau = mu1 - mu0
    
    # asymptotic variance of ATE
    phi1 = temp$A*(Ys-mu1)/temp$`piA`  # influence from rct treated
    phi0 = (1-temp$A)*(Ys-mu0)/(1-temp$`piA`)  # influence from rct control
    ## big Phi should be n by T
    Phi = cbind(phi1, phi0)
    B = (t(Phi) %*% Phi)/n
    ## sandwich
    sigma = B 
    
    # ATE as linear comb of parameters
    coef.mat = cbind(diag(T_follow), -diag(T_follow))
    # get standard error using the linear comb of var-cov matrix
    sd.tau = sqrt(diag(coef.mat %*% B %*% t(coef.mat)/(n)))
    
  }else {
    # propensity score model
    piS.model = glm(model_form_piS , data = df, family = "binomial")
    # outcome regression model
    Y0.model = lapply(model_form_mu0_ext, function(x){lm(as.formula(x), data = filter(df, A==0) )})
    Y0.model.dummy = lapply(model_form_mu0_ext, function(x){lm(as.formula(x), data = filter(df))})
    
    # predict Y0 from outcome regression models
    SuppressWarnings(Y0 = data.frame(sapply(1:T_follow, function(x){predict(Y0.model[[x]], newdata = df)})))
    colnames(Y0) = paste0("y", 1:T_follow, "_0")
    # for residual
    Yr = Y - Y0 
    colnames(Yr) = paste0("y", 1:T_follow, "_r")
    
    # estimate ATE
    SuppressWarnings(
    temp = df %>% cbind(., Y0, Yr) %>%
      mutate(piA = sum(A[S==1])/n,
             piS = sum(S)/(n+m),
             piSX = predict(piS.model, newdata = df, type="response"),
             rx = (piSX/(1-piSX))*((1-piS)/piS)) %>%
      mutate(w11 = piA, w10 = 1 - piA, w00 = rx)
    )
    ### create outcomes: obs * T
    Ys = as.matrix(Yr)
    
    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1 = colSums(potential[temp$S==1&temp$A==1,])/n
    mu10 = colSums(potential[temp$S==1&temp$A==0,])/n
    mu00 = colSums(potential[temp$S==0,])/(sum((1-temp$S)*temp$w00))
    
    
    ## Optimal weight as proposed in manuscript
    # fit1 = lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==1&A==0) )
    # sigma10 = (summary(fit1)$sigma)**2
    num = sum(temp$S*(1-temp$A)/temp$w10^2/(sum(temp$S*(1-temp$A)/temp$w10))^2)
    
    #+ sum(temp$S*(1-temp$A))*var(temp$S*(1-temp$A)*predict(fit1,newdata =temp)/temp$w10/sum(temp$S*(1-temp$A)/temp$w10))
    
    # fit0=lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==0) )
    # sigma00=(summary(fit0)$sigma)**2
    denom = sum((1-temp$S)*temp$w00^2/(sum((1-temp$S)*temp$w00))^2)
    #+ sum((1-temp$S))*var((1-temp$S)*predict(fit0,newdata =temp)*temp$w00/sum((1-temp$S)*temp$w00))
    w.opt = num/(num+denom)
    if(optimal_weight_flag){
      wt = w.opt
    }
    
    # final hybrid estimate as combination of RCT control and external control
    mu0 = (1-wt)*mu10 + wt*mu00
    tau = mu1 - mu0
    
    names(tau) = paste0("tau", 1:T_follow)
    
  }
  
  return(tau)
  
}