#' Using IPW with external borrowing
#'
#' @param model_form_piS 
#' @param optimal_weight_flag 
#' @param wt 
#' @param Bootstrap 
#' @param data 
#' @param outcome_col_name 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param covariates_col_name 
#' @param R 
#' @param bootstrap_CI_type 
#' @param alpha 
#' 
#' 
#' @include EC_IPW_OPT_bootstrap.R 
#' @return a list containing: tau (effect size), sd.tau (standard deviation), wt (weight)
#' @export
#'
#' @examples
#' EC_IPW_OPT(data = data,
#'    outcome_col_name = outcome_col_name, 
#'    trial_status_col_name = trial_status_col_name, 
#'    treatment_col_name = treatment_col_name, 
#'    covariates_col_name = covariates_col_name, 
#'    model_form_piS = model_form_piS,
#'    wt = wt, 
#'    optimal_weight_flag = optimal_weight_flag,
#'    Bootstrap, R, bootstrap_CI_type)

EC_IPW_OPT = function(data,
                      outcome_col_name, 
                      trial_status_col_name, 
                      treatment_col_name, 
                      covariates_col_name,
                      model_form_piS = "",
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
  names(S) = "S"
  names(A) = "A"
  model_form_piS = unlist(strsplit(model_form_piS, "~"))
  model_form_piS[1] = "S"
  model_form_piS = paste0(model_form_piS, collapse = "~")
  
  df = data.frame(Y, S = S, A = A, X)
  N = nrow(df)
  T_follow = ncol(Y)
  n = sum(df$S) # RCT sample size
  m = sum(1 - df$S) # external control sample size
  pi.S = n/N
  
  ## TODO: write running messages
  if (!quiet){
    cat("running IPW... \n")
  }

  
  if((!optimal_weight_flag) && wt == 0){

    # estimate ATE
    ## TODO: why do use the true propensity score? 
    temp = df %>%
      filter(S == 1) %>%
      mutate(piA = sum(A)/n) %>%
      mutate(w11 = piA, w10 = 1 - piA)
    
    ### create outcomes: obs by T
    Ys = as.matrix(Y[S==1, ])
    
    potential = data.frame((temp$A/temp$w11+(1-temp$A)/temp$w10)*Ys)
    mu1 = colSums(potential[temp$A==1, , drop = F])/n
    mu0 = colSums(potential[temp$A==0, , drop = F])/n
    
    # print(c(mu1, mu0))
    tau = mu1 - mu0
    names(tau) = paste0("tau", 1:T_follow)
    
    # asymptotic variance of ATE
    phi1 = temp$A*sweep(Ys, 2, mu1)/temp$`piA`  # influence from rct treated
    phi0 = (1-temp$A)*sweep(Ys, 2, mu0)/(1-temp$`piA`)  # influence from rct control
    ## big Phi should be n by T
    Phi = cbind(phi1, phi0)
    B = (t(Phi) %*% Phi)/n
    ## sandwich
    sigma = B 
    
    # ATE as linear comb of parameters
    coef.mat = cbind(diag(T_follow), -diag(T_follow))
    # get standard error using the linear comb of var-cov matrix
    sd.tau = sqrt(diag(coef.mat %*% B %*% t(coef.mat)/(n)))
    names(sd.tau) = paste0("sd.tau", 1:T_follow)
    
  }else{
    # propensity score model
    
    # TODO: generalize the following line to other possible models
    # piS.model = glm(as.formula(paste("S", "~", form_x)) , data = df, family = "binomial")
    piS.model = glm(as.formula(model_form_piS), data = df, family = "binomial")
    
    # print(piS.model)
    
    # estimate ATE
    temp = df%>%
      mutate(piA = sum(A[S==1])/n,
             piS = sum(S)/(n+m),
             piSX = predict(piS.model, newdata = df, type = "response"),
             # piSX = exp(log(2) + df$X - 3 * df$U)/(1+exp(log(2) + df$X - 3 * df$U)),
             rx = (piSX/(1 - piSX))*((1 - piS)/piS)) %>%
      mutate(w11 = piA, w10 = 1 - piA, w00 = rx)
    
    # temp$w00[temp$w00 > 0.9] = 0.9
    # temp$w00[temp$w00 < 0.1] = 0.1
    
    # print(max(abs(temp$w00)))
    
    ### create outcomes: obs * T
    Ys = as.matrix(Y)
    
    potential = data.frame((temp$S*temp$A/temp$w11 + temp$S*(1-temp$A)/temp$w10 + (1-temp$S)*temp$w00)*Ys)
    mu1 = colSums(potential[temp$S==1&temp$A==1, ,drop = F])/sum(temp$S*temp$A/temp$w11)
    mu10 = colSums(potential[temp$S==1&temp$A==0, ,drop = F])/sum(temp$S*(1-temp$A)/temp$w10)
    mu00 = colSums(potential[temp$S==0, ,drop = F])/(sum((1-temp$S)*temp$w00))
    
    
    # variance
    
    ## bread

    A11 = diag(rep(-1, T_follow), nrow = T_follow)
    A22 = diag(rep(-1, T_follow), nrow = T_follow)
    A33 = diag(rep(-mean((1-temp$S) * temp$w00/(1-temp$piS)), T_follow), nrow = T_follow)
    A34 = t((as.vector((1-temp$S) * temp$piSX/(temp$piS * (1-temp$piSX))) * sweep(Ys, 2, mu00))) %*% model.matrix(piS.model)/(n+m)
    piS.beta = t(model.matrix(piS.model)) %*% diag(-temp$`piSX`*(1-temp$`piSX`)) %*% model.matrix(piS.model)/(n+m)
    A44 = piS.beta
    
    n1=dim(A11)[1]; n2=dim(A22)[1]; n3=dim(A33)[1]; n4=dim(A44)[1]
    
    A = matrix(0,nrow=n1+n2+n3+n4,ncol=n1+n2+n3+n4)
    A[1:n1,1:n1]=A11
    A[(n1+1):(n1+n2),(n1+1):(n1+n2)]=A22
    A[(n1+n2+1):(n1+n2+n3),(n1+n2+1):(n1+n2+n3)]=A33
    A[(n1+n2+n3+1):(n1+n2+n3+n4),(n1+n2+n3+1):(n1+n2+n3+n4)]=A44
    
    A[(n1+n2+1):(n1+n2+n3),(n1+n2+n3+1):(n1+n2+n3+n4)]=A34
    
    
    ## meat
    phi1 = temp$S*temp$A*sweep(Ys, 2, mu1)/temp$piA/temp$piS    # influence from rct treated
    phi2 = temp$S*(1-temp$A)*sweep(Ys, 2, mu10)/(1-temp$piA)/temp$piS    # influence from rct control
    phi3 = (1-temp$S)*sweep(Ys, 2, mu00)*temp$w00/(1-temp$piS)
    phi.piS = (temp$S - temp$`piSX`)*model.matrix(piS.model)
    
    ## big Phi should be n by T+T+p,
    Phi = cbind(phi1, phi2, phi3, phi.piS)
    
    B = (t(Phi)%*%Phi)/(n+m)
    
    ## sandwich
    sigma = solve(A) %*% B %*% t(solve(A))
    
    ## Optimal weight as proposed in manuscript
    # fit1 = lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==1&A==0) )
    # sigma10 = (summary(fit1)$sigma)**2
    num = sum(temp$S*(1-temp$A)/temp$w10^2/(sum(temp$S*(1-temp$A)/temp$w10))^2)
    #+ sum(temp$S*(1-temp$A))*var(temp$S*(1-temp$A)*predict(fit1,newdata =temp)/temp$w10/sum(temp$S*(1-temp$A)/temp$w10))
    
    # fit0 = lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==0) )
    # sigma00 = (summary(fit0)$sigma)**2
    denom = sum((1-temp$S)*temp$w00^2/(sum((1-temp$S)*temp$w00))^2)
    #+ sum((1-temp$S))*var((1-temp$S)*predict(fit0,newdata =temp)*temp$w00/sum((1-temp$S)*temp$w00))
    w.opt = num/(num+denom)
    # only use the optimal weight if we want
    if(optimal_weight_flag == T){
      wt = w.opt
      # print(wt)
    }
    
    # final hybrid estimate as combination of RCT control and external control
    mu0 = (1-wt)*mu10 + wt*mu00
    tau = mu1-mu0
    
    names(tau) = paste0("tau", 1:T_follow)
    
    # ATE as linear comb of parameters
    coef.mat = cbind(diag(T_follow), 
                     -(1-wt)*diag(T_follow), 
                     -wt*diag(T_follow), 
                     matrix(0, nrow = T_follow, ncol = dim(piS.beta)[2]))
    # get standard error using the linear comb of var-cov matrix
    sd.tau = sqrt(diag(coef.mat%*% sigma %*%t(coef.mat)/(n+m)))
    
    names(sd.tau) = paste0("sd.tau", 1:T_follow)
    
  }
  
  # summarize results
  
  cutoff = qnorm(1-alpha/2, lower.tail = T)
  
  if(Bootstrap == T){
    Group_ID = df %>% group_by(S, A) %>% mutate(group_id = cur_group_id())
    Group_ID = Group_ID$group_id
    
    boot.ci.type = switch(bootstrap_CI_type,
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