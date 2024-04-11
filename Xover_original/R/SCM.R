#' Implement the Synthetic Control Method, in the manuscript  
#'
#' SCM() is the main function calculates the estimated ATE by SC method and Bootstrap CI, it calls subject_SC() and lambdacv().
#'
#' @param df A data frame contains simulated data, returned by cal.simXY()
#' @param method:
#' "DID-EC-OR" 
#' "DID-EC-IPW"
#' "DID-EC-AIPW"
#' @param setting The setting number
#' Setting 1: No unmeasured confounder 
#' Setting 2: unmeasured confounding with time-invariant effect, all models with observed covariates are correctly specified
#' Setting 3: unmeasured confounding with time-invariant effect, outcome models with observed covariates are correctly specified
#' Setting 4: unmeasured confounding with time-invariant effect, PS models with observed covariates are correctly specified
#' Setting 5: unmeasured confounding with time-invariant effect, both models with observed covariates are slightly misspecified
#' Setting 6: unmeasured confounding with time-variant effect, both models with observed covariates are slightly misspecified
#' @return A list contains: estimated ATE, SE, weihgt used, SE by Bootstrap and a 95% confidence interval for primary endpoint (only when Bootstrap=TRUE)
#' @examples
#' SCM(df=mysimXY)
#'  
SCM = function(df,
              setting=1){
  
  n0=sum(df$S==1 & df$A==0) # number of rct control subjects
  m=sum(df$S==0) # number of external control subjects
  
  # standardize numeric vars
  columns_to_standardize=c("SMN2_Copy_Number","Age_Enrollment","Y0")
  df[columns_to_standardize] <- lapply(df[columns_to_standardize], scale)
   
  # create data matrices: attributes by row and subject by column
  if(setting==1){
    X10 = t(as.matrix(df%>%filter(S==1&A==0)%>%dplyr::select(SMA_Type, SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    X00 = t(as.matrix(df%>%filter(S==0)%>%dplyr::select(SMA_Type, SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    
  }else{
    X10 = t(as.matrix(df%>%filter(S==1&A==0)%>%dplyr::select(SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    X00 = t(as.matrix(df%>%filter(S==0)%>%dplyr::select(SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    
  }
  
  # remove colnames of X10 and X00
  colnames(X00)=NULL
  colnames(X10)=NULL
  
  # Use LOOCV to find the optimal penalty parameter lambda
  lambda=lambdacv(setting=1,ec=X00,lower=-60,upper=0.000001,trials=10)
  
  # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
  res = lapply(1:n0, subject_SC, setting=setting, X10=X10,X00=X00,lambda=lambda)

  # weight matrix: n0*m (combinging by rows)
  wt.mat=do.call(rbind,lapply(res, function(x) as.vector(x[[1]])))

  # estimated outcome matrix: n0*T (combining by rows)
  y.est.mat=do.call(rbind,lapply(res, function(x) as.vector(x[[2]])))

  # Aggregate group level synthetic control estiamte
  Y.trt=filter(df, S==1&A==1)%>%
    summarise(Y3obs=mean(Y3),Y4obs=mean(Y4))%>%
    dplyr::select(Y3obs, Y4obs)%>%unlist()
  tau=Y.trt-colMeans(y.est.mat)

  ####### Use Bootstrap for standard error and confidence intervals

  boot.out <- boot(data=df,
                   statistic=SCMboot,
                   R=100,
                   # parallel = "multicore",
                   # ncpus=4,
                   setting=setting,
                   lambda=lambda)

  return(list(c(tau[1], sd(boot.out$t[,1]), boot.ci(boot.out,index=1,type=c("perc"))$percent[4:5]),
              c(tau[2], sd(boot.out$t[,2]), boot.ci(boot.out,index=2,type=c("perc"))$percent[4:5])))
  
}


#' find synthetic control for one specific subject

subject_SC =function(setting,subject,X10,X00,lambda){
  if(setting==1){
    match.vars=which(rownames(X10) %in% c("SMA_Type", "SMN2_Copy_Number", "Scoliosis", "Age_Enrollment", "Y0", "Y1", "Y2"))
  }else{
    match.vars=which(rownames(X10) %in% c("SMN2_Copy_Number", "Scoliosis", "Age_Enrollment", "Y0", "Y1", "Y2"))
  }
  # subset to only the intersted subject and matching vars
  x1 = X10[match.vars,subject]
  X0 = X00[match.vars,]
  
  # optimization
  w<-Variable(dim(X00)[2])
  loss <- sum(((x1 - X0%*%w))^2)
  penal<-function(w, x1, X0, lambda=lambda){
    lambda*(sum(colSums((x1-X0)^2)*w))
    }
  obj<-loss+penal(w=w,x1=x1,X0=X0,lambda=lambda)
  constr <- list(sum(w) == 1, w>= 0)
  prob <- Problem(Minimize(obj), constr)
  result <- solve(prob, solver='ECOS') # OSQP, SCS, ECOS
  # estimated weight for this subject
  wt.est=result$getValue(w)
  # sc estimate for all time points for this subject
  y.est=X00[c("Y3","Y4"),]%*%wt.est
  
  return(list(wt.est,y.est))
}

#' find the optima lambda via LOOCV 

lambdacv = function(setting=1,
                   ec,
                   lower=-20,
                   upper=0.1,
                   trials=10){
  
 
  lambda_vals = 10^seq(lower, log10(upper), length.out = trials)
  
  # leave-one-out cv error for individul external controls
  mse_vals = sapply(lambda_vals,
                     function (lambda) {
                       res=lapply(1:dim(ec)[2], function(loocv){
                          
                         if(setting==1){
                           match.vars=which(rownames(ec) %in% c("SMA_Type", "SMN2_Copy_Number", "Scoliosis", "Age_Enrollment", "Y0", "Y1", "Y2"))
                         }else{
                           match.vars=which(rownames(ec) %in% c("SMN2_Copy_Number", "Scoliosis", "Age_Enrollment", "Y0", "Y1", "Y2"))
                         }
                         # subset to only the intersted subject and matching vars
                         x1 = ec[match.vars,loocv]
                         X0 = ec[match.vars,-loocv]
                         
                         # optimization
                         w<-Variable(dim(ec)[2]-1)
                         loss <- sum(((x1 - X0%*%w))^2)
                         penal<-function(w,x1,X0,lambda=lambda){lambda*(sum(colSums((x1-X0)^2)*w))}
                         obj<-loss+penal(w=w,x1=x1,X0=X0,lambda=lambda)
                         constr <- list(sum(w) == 1, w>= 0)
                         prob <- Problem(Minimize(obj),constr)
                         result <- solve(prob,solver='ECOS') # OSQP, SCS, ECOS
                         # estimated weight for this subject
                         wt.est=result$getValue(w)
                         # sc estimate for all time points for this subject
                         y.est=ec[c("Y3","Y4"),-loocv]%*%wt.est
                         
                         return(list(wt.est,y.est))
                       })
                       # estimated outcome matrix: n0*T (combining by rows)
                       y.est.mat=do.call(rbind,lapply(res, function(x) as.vector(x[[2]])))
                       
                       # calcualte MSE
                       mean((ec[c("Y3","Y4"),]-t(y.est.mat))^2)
                       
                     })
  
  return(lambda_vals[which.min(mse_vals)])
  
}
 
