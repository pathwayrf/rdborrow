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
#' system.time(res<-SCM(df=mysimXY))
#'  
SCM<-function(df=mysimXY,
              setting=1){
  
  # change X3 to be numeric for matching
  df=df%>%mutate(X1=ifelse(X1==1,1,0),X2=ifelse(X2==1,1,0))
  
  n0=sum(df$S==1 & df$A==0) # number of rct control subjects
  m=sum(df$S==0) # number of external control subjects
  
  # standardize numeric vars
  if(setting %in% c(1,2,3,4)){
    columns_to_standardize=c("X3","X4","X5","Y1","Y2")
  }else{
    columns_to_standardize=c("X3","X4","X5","Y1","Y2","Y3","Y4")
  }
  df[columns_to_standardize] <- lapply(df[columns_to_standardize], scale)

  
  # create data matrices: attributes by row and subject by column
  if(setting %in% c(1,2,3,4)){
    X10 = t(as.matrix(df%>%filter(S==1&A==0)%>%dplyr::select(X1:X5,Y1:Y4)))
    X00 = t(as.matrix(df%>%filter(S==0)%>%dplyr::select(X1:X5,Y1:Y4)))
    
  }else{
    X10 = t(as.matrix(df%>%filter(S==1&A==0)%>%dplyr::select(X1:X5,Y1:Y6)))
    X00 = t(as.matrix(df%>%filter(S==0)%>%dplyr::select(X1:X5,Y1:Y6)))
    
  }
  
  
  # remove colnames of X10 and X00
  colnames(X00)=NULL
  colnames(X10)=NULL

  
  # Use LOOCV to find the optimal penalty parameter lambda
  lambda=lambdacv(setting=1,ec=X00,lower=-10,upper=1,trials=10)
  
  # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
  #n_cores <- detectCores()
  res=lapply(1:n0, subject_SC, setting=setting, X10=X10,X00=X00,lambda=lambda)
   
  # weight matrix: n0*m (combinging by rows)
  wt.mat=do.call(rbind,lapply(res, function(x) as.vector(x[[1]])))

  # estimated outcome matrix: n0*T (combining by rows)
  y.est.mat=do.call(rbind,lapply(res, function(x) as.vector(x[[2]])))

  # Aggregate group level synthetic control estiamte
  if(setting %in% c(1,2,3,4)){
    Y.trt=filter(df, S==1&A==1)%>%
      summarise(Y3obs=mean(Y3),Y4obs=mean(Y4))%>%
      dplyr::select(Y3obs, Y4obs)%>%unlist()
    tau=Y.trt-colMeans(y.est.mat)
  }else{
    Y.trt=filter(df, S==1&A==1)%>%
      summarise(Y5obs=mean(Y5),Y6obs=mean(Y6))%>%
      dplyr::select(Y5obs, Y6obs)%>%unlist()
    tau=Y.trt-colMeans(y.est.mat)
  }

  ####### Use Bootstrap for standard error and confidence intervals

  # boot.out <- boot(data=df,
  #                  statistic=SCMboot,
  #                  R=5,
  #                  strata = df$S,
  #                  parallel = "multicore",
  #                  # ncpus=4,
  #                  setting=setting,
  #                  lambda=lambda)
  # 
  # return(list(c(tau[1], sd(boot.out$t[,1]), boot.ci(boot.out,index=1,type=c("perc"))$percent[4:5]),
  #             c(tau[2], sd(boot.out$t[,2]), boot.ci(boot.out,index=2,type=c("perc"))$percent[4:5])))
  # 
  return(list(c(tau[1], rep(NA,3)),
              c(tau[2], rep(NA,3))))
  

  # # random sample index with replacement
  # ind.list=replicate(n = 100, 
  #                    expr = {sample(1:nrow(df), size=nrow(df), replace = TRUE)},simplify =FALSE)
  # 
  # n_cores <- detectCores()
  # ans.list=mclapply(ind.list, function(index){
  #   tryCatch(
  #     expr = SCMboot(index),
  #     error = function(e) return(NULL) # return NULL or some other value on error
  #   )
  # }, mc.cores = n_cores)
  # 
  # tmp=do.call(rbind, ans.list)
  # cis=apply(tmp,2,quantile, probs = c(0.025, 0.975))
  # ses=apply(tmp,2,sd)

  
  # return(list(c(tau[1], ses[1], cis[,1]),
  #             c(tau[2], ses[2], cis[,2])))


}

# SCMboot<-function(index){
#   df_b=df[index,]
#   
#   n0=sum(df_b$S==1 & df_b$A==0) # number of rct control subjects
#   m=sum(df_b$S==0) # number of external control subjects
#   
#   
#   # create data matrices: attributes by row and subject by column
#   if(setting %in% c(1,2,3,4)){
#     X10 = t(as.matrix(df_b%>%filter(S==1&A==0)%>%dplyr::select(X1:X5,Y1:Y4)))
#     X00 = t(as.matrix(df_b%>%filter(S==0)%>%dplyr::select(X1:X5,Y1:Y4)))
#     
#   }else{
#     X10 = t(as.matrix(df_b%>%filter(S==1&A==0)%>%dplyr::select(X1:X5,Y1:Y6)))
#     X00 = t(as.matrix(df_b%>%filter(S==0)%>%dplyr::select(X1:X5,Y1:Y6)))
#     
#   }
#   
#   # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
#   #n_cores <- detectCores()
#   ans=lapply(1:n0, subject_SC, setting=setting, X10=X10,X00=X00,lambda=lambda)
#   
#   # estimated outcome matrix: n0*T (combining by rows)
#   y.est.mat=do.call(rbind,lapply(ans, function(x) as.vector(x[[2]])))
#   
#   # Aggregate group level synthetic control estiamte
#   if(setting %in% c(1,2,3,4)){
#     Y.trt=filter(df_b, S==1&A==1)%>%
#       summarise(Y3obs=mean(Y3),Y4obs=mean(Y4))%>%
#       dplyr::select(Y3obs, Y4obs)%>%unlist()
#     tau=Y.trt-colMeans(y.est.mat)
#   }else{
#     Y.trt=filter(df_b, S==1&A==1)%>%
#       summarise(Y5obs=mean(Y5),Y6obs=mean(Y6))%>%
#       dplyr::select(Y5obs, Y6obs)%>%unlist()
#     tau=Y.trt-colMeans(y.est.mat)
#   }
#   
#   return(tau)
# }
  

#' find synthetic control for one specific subject

subject_SC =function(setting,subject,X10,X00,lambda){
  if(setting %in% c(1)){
    match.vars=which(rownames(X00) %in% c("X1", "X2", "X3", "X4", "X5", "Y1", "Y2"))
  }else if(setting %in% c(2,3,4)){
    match.vars=which(rownames(X00) %in% c("X2", "X3", "X4", "X5", "Y1", "Y2"))
  }else{
    match.vars=which(rownames(X00) %in% c("X2", "X3", "X4", "X5", "Y1", "Y2", "Y3", "Y4"))
  }
  # subset to only the intersted subject and matching vars
  x1 = X10[match.vars,subject]
  X0 = X00[match.vars,]
  
  # optimization
  w<-Variable(dim(X00)[2])
  loss <- sum(((x1 - X0%*%w))^2)
  penal<-function(w,x1,X0,lambda=lambda){lambda*(sum(colSums((x1-X0)^2)*w))}
  obj<-loss+penal(w=w,x1=x1,X0=X0,lambda=lambda)
  constr <- list(sum(w) == 1, w>= 0)
  prob <- Problem(Minimize(obj),constr)
  result <- solve(prob,solver='ECOS') # OSQP, SCS, ECOS
  # estimated weight for this subject
  wt.est=result$getValue(w)
  # sc estimate for all time points for this subject
  if(setting %in% c(1,2,3,4)){
    y.est=X00[c("Y3","Y4"),]%*%wt.est
     
  }else{
    y.est=X00[c("Y5","Y6"), ]%*%wt.est
  }
  
   
  
  return(list(wt.est,y.est))
}

#' find the optima lambda via LOOCV 

lambdacv<-function(setting=1,
                   ec,
                   lower=-10,
                   upper=1,
                   trials=10){
  
 
  lambda_vals <- 10^seq(lower, log10(upper), length.out = trials)
  
  # Detect the number of cores in your machine
  n_cores <- detectCores()
  
  # leave-one-out cv error for individul external controls
  mse_vals <- mclapply(lambda_vals,
                      function (lambda) {
                        res=lapply(1:dim(ec)[2], function(loocv){
                          
                          if(setting %in% c(1)){
                            match.vars=which(rownames(ec) %in% c("X1", "X2", "X3", "X4", "X5", "Y1", "Y2"))
                          }else if(setting %in% c(2,3,4)){
                            match.vars=which(rownames(ec) %in% c( "X2", "X3", "X4", "X5", "Y1", "Y2"))
                          }else{
                            match.vars=which(rownames(ec) %in% c( "X2", "X3", "X4", "X5", "Y1", "Y2", "Y3", "Y4"))
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
                          if(setting %in% c(1,2,3,4)){
                            y.est=ec[c("Y3","Y4"),-loocv]%*%wt.est
                          }else{
                            y.est=ec[c("Y5","Y6"),-loocv]%*%wt.est
                          }
                          
                          
                          return(list(wt.est,y.est))
                        })
                        # estimated outcome matrix: n0*T (combining by rows)
                        y.est.mat=do.call(rbind,lapply(res, function(x) as.vector(x[[2]])))
                        
                        # calcualte MSE
                        if(setting %in% c(1,2,3,4)){
                          mean((ec[c("Y3","Y4"),]-t(y.est.mat))^2)
                        }else{
                          mean((ec[c("Y5","Y6"),]-t(y.est.mat))^2)
                        }
                        
                        
                      }, mc.cores = n_cores)
  
   
  return(lambda_vals[which.min(do.call(c,mse_vals))])
  
}
 
