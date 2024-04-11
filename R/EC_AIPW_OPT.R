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
#' @param Bootstrap 
#' @param R 
#' @param bootstrap_CI_type 
#' @param alpha 
#'
#' @include EC_AIPW_OPT_bootstrap.R
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

EC_AIPW_OPT = function(data,
                       outcome_col_name, 
                       trial_status_col_name, 
                       treatment_col_name, 
                       covariates_col_name,
                       model_form_piS = "",
                       model_form_mu0_ext = "",
                       optimal_weight_flag = F,
                       wt = 0,
                       Bootstrap = F,
                       R = 5e2,
                       bootstrap_CI_type = "bca",
                       alpha = 0.05,
                       quiet = TRUE){
  Y = subset(data, select = outcome_col_name)
  S = subset(data, select = trial_status_col_name)
  A = subset(data, select = treatment_col_name)
  X = subset(data, select = covariates_col_name)
  
  df = data
  N = nrow(df)
  T_follow = ncol(Y)
  n = sum(df$S) # RCT sample size
  m = sum(1 - df$S) # external control sample size
  pi.S = n/N
  
  ## TODO: write running messages
  if(!quiet){
    cat("running AIPW... \n")
  }
  
  # if((!optimal_weight_flag) && wt == 0){
  #   # estimate ATE
  #   ## TODO: why do use the true propensity score? 
  #   temp = df %>%
  #     filter(S==1) %>%
  #     mutate(`piA`=sum(A)/n) %>%
  #     mutate(w11 = `piA`, w10 = 1 - `piA`)
  #   
  #   ### create outcomes: obs by T
  #   Ys = as.matrix(Y[S==1, ])
  #   
  #   potential = (temp$A/temp$w11+(1-temp$A)/temp$w10)*Ys
  #   mu1 = colSums(potential[temp$A==1, ])/n
  #   mu0 = colSums(potential[temp$A==0, ])/n
  #   
  #   tau = mu1 - mu0
  #   
  #   # asymptotic variance of ATE
  #   phi1 = temp$A*(Ys-mu1)/temp$`piA`  # influence from rct treated
  #   phi0 = (1-temp$A)*(Ys-mu0)/(1-temp$`piA`)  # influence from rct control
  #   ## big Phi should be n by T
  #   Phi = cbind(phi1, phi0)
  #   B = (t(Phi) %*% Phi)/n
  #   ## sandwich
  #   sigma = B 
  #   
  #   # ATE as linear comb of parameters
  #   coef.mat = cbind(diag(T_follow), -diag(T_follow))
  #   # get standard error using the linear comb of var-cov matrix
  #   sd.tau = sqrt(diag(coef.mat %*% B %*% t(coef.mat)/(n)))
  #   
  # }else 
    {
    # propensity score model
    piS.model = glm(model_form_piS , data = df, family = "binomial")
    # outcome regression model
    Y0.model = lapply(model_form_mu0_ext, function(x){lm(as.formula(x), data = filter(df, A==0) )})
    Y0.model.dummy = lapply(model_form_mu0_ext, function(x){lm(as.formula(x), data = filter(df))})
    
    # predict Y0 from outcome regression models
    Y0 = data.frame(sapply(1:T_follow, function(x){predict(Y0.model[[x]], newdata = df)}))
    colnames(Y0) = paste0("y", 1:T_follow, "_0")
    # for residual
    Yr = Y - Y0 
    colnames(Yr) = paste0("y", 1:T_follow, "_r")
    
    # estimate ATE
    temp = df %>% cbind(., Y0, Yr) %>%
      mutate(piA = sum(A[S==1])/n,
             piS = sum(S)/(n+m),
             piSX = predict(piS.model, newdata = df, type="response"),
             rx = (piSX/(1-piSX))*((1-piS)/piS)) %>%
      mutate(w11 = piA, w10 = 1 - piA, w00 = rx)
    
    ### create outcomes: obs * T
    Ys = as.matrix(Yr)
    
    
    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1 = colSums(potential[temp$S==1&temp$A==1, , drop = F])/n
    mu10 = colSums(potential[temp$S==1&temp$A==0, , drop = F])/n
    mu00 = colSums(potential[temp$S==0, , drop = F])/(sum((1-temp$S)*temp$w00))
    
    
    
    # variance
    
    ## bread
    A11 = diag(rep(-1, T_follow), nrow = T_follow)
    A22 = diag(rep(-1, T_follow), nrow = T_follow)
    A33 = diag(rep(-mean((1-temp$S)*temp$w00/(1-temp$piS)), T_follow), nrow = T_follow)
    A34 = t((as.vector((1-temp$S)*temp$`piSX`/(temp$`piS`*(1-temp$`piSX`)))*(Ys-mu00)))%*%model.matrix(piS.model)/(n+m)
    piS.beta = t(model.matrix(piS.model))%*%diag(-temp$`piSX`*(1-temp$`piSX`))%*%model.matrix(piS.model)/(n+m)
    A44 = piS.beta
    
    n1 = dim(A11)[1]
    n2 = dim(A22)[1]
    n3 = dim(A33)[1]
    n4 = dim(A44)[1]
    
    A0 = as.matrix(bdiag(A11, A22, A33, A44))
    A0[(n1+n2+1):(n1+n2+n3),(n1+n2+n3+1):(n1+n2+n3+n4)] = A34
    
    
    Y0.model.dim = sapply(1:T_follow, function(x){dim(model.matrix(Y0.model.dummy[[x]]))[2]})
    n5 = sum(Y0.model.dim)
    
    A45 = matrix(0, nrow = n4, ncol = n5)
    A51 = matrix(0, nrow = n5, ncol = n1+n2+n3+n4)
    
    # cat(length(as.vector(-temp$S*temp$A/(temp$piS*(temp$piA)))))
    # cat(dim(model.matrix(Y0.model.dummy[[1]])))
    
    Phi1.gamma.list = lapply(1:T_follow, 
                             function(x){
                               as.vector(-temp$S*temp$A/(temp$piS*(temp$piA))) %*% model.matrix(Y0.model.dummy[[x]])/(n+m)
                             })
    Phi1.gamma = as.matrix(do.call(bdiag, Phi1.gamma.list))
    
    Phi2.gamma.list = lapply(1:T_follow, 
                             function(x){
                               as.vector(-temp$S*(1-temp$A)/(temp$piS*(1-temp$piA))) %*% model.matrix(Y0.model.dummy[[x]])/(n+m)
                             })
    Phi2.gamma = as.matrix(do.call(bdiag, Phi2.gamma.list))
    
    Phi3.gamma.list = lapply(1:T_follow, 
                             function(x){
                               as.vector(-(1-temp$S)*temp$rx/(1-temp$piS)) %*% model.matrix(Y0.model.dummy[[x]])/(n+m)
                             })
    Phi3.gamma = as.matrix(do.call(bdiag, Phi3.gamma.list))
    
    Y0.gamma.list = lapply(1:T_follow, 
                           function(x){ 
                             -t(model.matrix(Y0.model.dummy[[x]])) %*% diag((1-temp$A)/(1-mean(temp$A))) %*% model.matrix(Y0.model.dummy[[x]])/(n+m)
                           })
    Y0.gamma = as.matrix(do.call(bdiag, Y0.gamma.list))
    
    
    # cat(c(n1, n2, n3, n4, n5))
    A.left = rbind(A0, A51)
    A.right = rbind(Phi1.gamma, Phi2.gamma, Phi3.gamma, A45, Y0.gamma)
    A = cbind(A.left, A.right)
    
    ## meat
    phi1 = temp$S*temp$A*(Ys-mu1)/temp$`piA`/temp$piS    # influence from rct treated
    phi2 = temp$S*(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)/temp$piS    # influence from rct control
    phi3 = (1-temp$S)*(Ys-mu00)*temp$w00/(1-temp$piS)
    phi.piS = (temp$S-temp$`piSX`)*model.matrix(piS.model)
    phi.Y0 = do.call(cbind, 
                     args = lapply(1:T_follow, function(x){
                       ((1-temp$A)/(1-mean(temp$A))) * (Ys[, x] * model.matrix(Y0.model.dummy[[x]]))
                     }))
    
    ## big Phi should be n by 4+2+p, observations of same phi in the same column
    Phi = cbind(phi1, phi2, phi3, phi.piS, phi.Y0)
    
    B = (t(Phi)%*%Phi)/(n+m)
    
    
    ## sandwich
    sigma = solve(A)%*%B%*%t(solve(A))
    
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
    if(optimal_weight_flag){wt = w.opt}
    
    # final hybrid estimate as combination of RCT control and external control
    mu0 = (1-wt)*mu10 + wt*mu00
    tau = mu1 - mu0
    
    names(tau) = paste0("tau", 1:T_follow)
    
    # ATE as linear comb of parameters
    ## TODO: too many zeros; maybe can be improved?
    coef.mat = cbind(diag(T_follow), 
                     -(1-wt)*diag(T_follow), 
                     -wt*diag(T_follow), 
                     matrix(0, nrow = T_follow, ncol = dim(piS.beta)[2] + n5))
    # get standard error using the linear comb of var-cov matrix
    sd.tau = sqrt(diag(coef.mat %*% sigma %*% t(coef.mat)/(n+m)))
    
    names(sd.tau) = paste0("sd.tau", 1:T_follow)
    
  }
  
  cutoff = qnorm(1-alpha/2, lower.tail = T)
  
  if(Bootstrap == T){
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
                     statistic = EC_IPW_OPT_bootstrap, 
                     outcome_col_name = outcome_col_name, 
                     trial_status_col_name = trial_status_col_name, 
                     treatment_col_name = treatment_col_name, 
                     covariates_col_name = covariates_col_name,
                     model_form_piS = model_form_piS,
                     optimal_weight_flag = optimal_weight_flag,
                     wt = wt,
                     R = R, 
                     strata = Group_ID)
    lower_CI_boot = sapply(1:T_follow, 
                           function(x) {
                             ci = boot.ci(boot.out, conf = 1 - alpha, type = bootstrap_CI_type, index = x)
                             ci[[boot.ci.type]][4]
                           })
    upper_CI_boot = sapply(1:T_follow, 
                           function(x) {
                             ci = boot.ci(boot.out, conf = 1 - alpha, type = bootstrap_CI_type, index = x)
                             ci[[boot.ci.type]][5]
                           })
    
    sd.boot=sqrt(diag(var(boot.out$t)))
    
    results = data.frame(
      point_estimates = tau,
      standard_deviation = sd.tau,
#      lower_CI_normal = tau - sd.tau * cutoff,
#      upper_CI_normal = tau + sd.tau * cutoff,
      lower_CI_boot = lower_CI_boot,
      upper_CI_boot = upper_CI_boot
    )
    return(list(results = results, borrow_weight = wt))
    
  }else{
    results = data.frame(
      point_estimates = tau,
      standard_deviation = sd.tau,
      lower_CI_normal = tau - sd.tau * cutoff,
      upper_CI_normal = tau + sd.tau * cutoff
    )
    return(list(results = results, borrow_weight = wt))
    
  }
  
}