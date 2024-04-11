#' Calibrate simulation according to SUNFISH and external data
#'
#' This function simulate X(baseline covariates), S(trial participation), A(treatment assignemnt), Y(repeated outcomes), according to the real data, and return a data frame and true ATE.
#'
#' @param path Path to where the proceed SUNFISH+external data, i.e. "bene.RData" and "outcome.RData", are located.
#' @param setting The setting number
#' Setting 1: all confounder correctly specified
#' Setting 2: outcome model misspecifed
#' Setting 3: propensity model misspecifed
#' Setting 4: both propensity model and outcome models misspecified
#' Setting 5: Exists unmeasured confounder: SMM2_Copy number as pseudo omitted confounder
#' @param size Sample size
#' @param form_x_s Specify the assumed true functional relation between baseline covariates and the propensity of trial parcipation
#' @param form_x Specify the assumed true functional relation between baseline covariates and the outcomes (at each time point)
#' @param form_x_s_mis For settings require mis-specified models, we use a slight departure from linear terms, such as quadratic terms, for propensity score model
#' @param form_x_mis For settings require mis-specified models, we use a slight departure from linear terms, such as quadratic terms, for outcome models
#' @return a list contains simulated data and true ATE
#' @examples
#' result=cal.simXY(setting=1)
#' list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)



cal.simXY=function(path=mydir,
                   setting=1,
                   size=1000){
  
  # load the pre-proceed data
  load(paste(path,"myproj2/output/bene.RData",sep = "")) # patient info
  load(paste(path,"myproj2/output/outcome.RData",sep = "")) # outcome measures: each row is patient-measure, so one patient has multiple lines corresponding to different measures at different times
  
  # merge the patient level file with the patient-outcome file, only keep the relevant outcome measures
  # the resulting data frame mydat should contain 221 patients in total, with 164 SUNFISH patients
  # a quik tabulation: table(mydat$S)[1]/sum(table(mydat$S)), shows: the internal treated: internal control: external control are roughly in 100:50:50 ratio
  mydat=bene%>%
    filter(Study %in% c("Previous Placebo","Sunfish Part 2"))%>%
    mutate(Ambulatory_Status=ifelse(Ambulatory_Status=="NON-AMBULATORY",0,ifelse(Ambulatory_Status=="AMBULATORY",1,NA)))%>%
    mutate(SMA_Type=ifelse(SMA_Type=="TYPE II",0,ifelse(SMA_Type=="TYPE III",1,NA)))%>%
    mutate(SMN2_Copy_Number=as.numeric(SMN2_Copy_Number))%>%
    mutate(Scoliosis=ifelse(Scoliosis=="Y",1,ifelse(Scoliosis=="N",0,NA)))%>%
    mutate(Age_Enrollment=Age)%>%
    mutate(S=ifelse(Study=="Sunfish Part 2",1,0))%>%
    mutate(A=ifelse(ACTARM=="PART 2 - RO7034067",1,0))%>%
    dplyr::select(UNI_ID,Study,S,A,Ambulatory_Status, SMA_Type, SMN2_Copy_Number, Scoliosis, Age_Enrollment)%>%
    left_join((outcome%>%filter(AVISIT=="BASELINE")%>%dplyr::select(UNI_ID,MFM32)%>%rename(Y0=MFM32)),by="UNI_ID")%>%
    left_join((outcome%>%filter(AVISIT=="Week 35"|AVISIT=="Week 26")%>%dplyr::select(UNI_ID,change_from_baseline)%>%rename(Y1=change_from_baseline)),by="UNI_ID")%>%
    left_join((outcome%>%filter(AVISIT=="Week 52")%>%dplyr::select(UNI_ID,change_from_baseline)%>%rename(Y2=change_from_baseline)),by="UNI_ID")%>%
    left_join((outcome%>%filter(AVISIT=="Week 78")%>%dplyr::select(UNI_ID,change_from_baseline)%>%rename(Y3=change_from_baseline)),by="UNI_ID")%>%
    left_join((outcome%>%filter(AVISIT=="Week 104")%>%dplyr::select(UNI_ID,change_from_baseline)%>%rename(Y4=change_from_baseline)),by="UNI_ID")
  
  
  # Generate baseline covariates by conditioning strategy, one by one
  # Order of conditioning: SMA_Type (0 or 1) -> Scoliosis (0 or 1) -> SMN2_Copy_Number (2 3 or 4) -> Age_Enrollment give all previous -> Y0  given all previous
  fit=lm(SMA_Type~1, data=mydat) 
  SMA_Type=sample (c(0,1), size=size, replace=T, prob=c(1-fit$coefficients,fit$coefficients))
  simXY=data.frame(id=1:size,SMA_Type)
  
  mylogit <- glm( Scoliosis ~SMA_Type, data = mydat, family = "binomial")
  pred <- predict(mylogit, newdata = simXY,type='response'  )
  simXY$Scoliosis=sapply(pred, rbinom, n=1,size=1)
   
  mylogit <- multinom(  SMN2_Copy_Number~SMA_Type+Scoliosis, data = mydat)
  pred <- predict(mylogit, newdata = simXY, type="probs")
  simXY$SMN2_Copy_Number=sapply(split(pred,row(pred)), function(prob)sample(x=c(2,3,4), size=1, prob=prob))
  
  #fit=lm(log(Age_Enrollment)~factor(SMN2_Copy_Number)+factor(SMA_Type)+factor(Scoliosis), data=mydat)
  # logmu <- predict(fit, newdata = simXY)
  # plot(density(log(mydat$Age_Enrollment)))
  # plot(density(logmu))
  #simXY$Age_Enrollment=exp(sapply(logmu, rnorm, n=1, sd=summary(fit)$sigma))
  simXY$Age_Enrollment=rnorm(n=size, mean=mean(mydat$Age_Enrollment), sd=sd(mydat$Age_Enrollment))
  # plot(density(simXY$Age_Enrollment))
  # plot(density(mydat$Age_Enrollment))


  #fit=lm(Y0~Age_Enrollment+SMN2_Copy_Number+SMA_Type+Scoliosis, data=mydat)
  fit=lm(Y0~SMN2_Copy_Number+SMA_Type+Scoliosis, data=mydat)
  mu <- predict(fit, newdata = simXY)
  simXY$Y0=sapply(mu, rnorm,n=1,sd=summary(fit)$sigma)
  
  
  ## rename SMA_Type (0 or 1) as X1-> Scoliosis (0 or 1) as X2 -> SMN2_Copy_Number (2 3 or 4) as X3 -> Age_Enrollment as X4 -> Y0 as X5
  simXY=simXY%>%rename(X1=SMA_Type,
                       X2=Scoliosis,
                       X3=SMN2_Copy_Number,
                       X4=Age_Enrollment,
                       X5=Y0)%>%mutate(X1=factor(X1),X2=factor(X2))

  mydat=mydat%>%rename(X1=SMA_Type,
                       X2=Scoliosis,
                       X3=SMN2_Copy_Number,
                       X4=Age_Enrollment,
                       X5=Y0)%>%mutate(X1=factor(X1),X2=factor(X2))

  # wide to long format
  lmydat=mydat%>%
    gather(key = "time", value="Y",Y1:Y4)%>%
    mutate(time=as.numeric(substr(time,2,2)))%>%
    arrange(UNI_ID, time)
  
  
  
  if(setting %in% 1:3){
    # 1: no unmeasured confounding, no time-varying coef Us~Ys necessary
    # 2: has unmeasured confounding U, no time-varying coef Us~Ys 
    
    # true ps model with pars chosen by data
    ps.model <- glm(S~X1+X2+X3+X4+X2*X4+X5, data = mydat, family = "binomial")
    out.model <- lm(Y~X1+X2+X3+X4+X2*X4+X5+A*time, data=lmydat)
    
    # change true ATE to be same across settings
    out.model$coefficients[c("A","A:time")]=c(0,2.5/4)
    
    # stronger spurious correlation
    out.model$coefficients["X21:X4"]=0.4
    #ps.model$coefficients["X21:X4"]=1.4
    
  }else if(setting %in% c(4)){
    # 3: has unmeasured confounding U and X1 as U, time-varying coef X1~Ys
     
    # true ps model with pars chosen by data
    ps.model <- glm(S~X1+X2+X3+X4+X2*X4+X5, data = mydat, family = "binomial")
    out.model <- lm(Y~X1+X2+X3+X4+X2*X4+X5+A*time+X1*time, data=lmydat)
    
    # change true ATE to be same across settings
    out.model$coefficients[c("A","A:time")]=c(0,2.5/4)
    
    # stronger spurious correlation
    out.model$coefficients["X21:X4"]=0.4
    #ps.model$coefficients["X21:X4"]=1.4
    
  } else if(setting %in% c(5)){
     
    # 4: same as 3 but let T1 consists of 4+ time points
    
    # true ps model with pars chosen by data
    ps.model <- glm(S~X1+X2+X3+X4+X2*X4+X5, data = mydat, family = "binomial")
    out.model <- lm(Y~X1+X2+X3+X4+X2*X4+X5+A*time+X1*time, data=lmydat)
    
    
    # change true ATE to be same across settings
    out.model$coefficients[c("A","A:time")]=c(0,2.5/6)
    
    # stronger spurious correlation
    out.model$coefficients["X21:X4"]=0.4
    #ps.model$coefficients["X21:X4"]=1.4
    
  }
   
  # Set residual variation  
  sigma_err=0.5 #summary(out.model)$sigma
  #sigma_b=0.1
  
   
  # simulate S
  pred <-predict(ps.model, newdata =simXY,type="response")
  S=sapply(pred, rbinom, n=1, size=1)
  
  # simulate A: among S==1 randomization!
  piA=mean(filter(mydat,S==1)$A)
  A=S*rbinom(n=size,size=1,prob=piA)
  study=ifelse(S==1 & A==0,"RCT Control Arm",ifelse(S==1 & A==1,"RCT Trt Arm","External Control"))
  simXY=cbind(simXY,study,S,A) 
  
   
  
  # Simulate outcome Ys
  if(setting %in% c(1,2,4)){
    lsimXY=simXY%>%
      left_join(data.frame(id=rep(simXY$id, each=4), time=rep(1:4, size)), by="id")
    
    lsimXY=lsimXY%>%
      mutate(Y=sapply(predict(out.model, newdata = lsimXY), rnorm, n=1, sd=sigma_err))
    
  }else if(setting==3){
    lsimXY=simXY%>%
      left_join(data.frame(id=rep(simXY$id, each=4), time=rep(1:4, size)), by="id")
    
    # study bias set to 0.25 =50% of incremetanl ATE
    lsimXY=lsimXY%>%
      mutate(Y=sapply(predict(out.model, newdata = lsimXY), rnorm, n=1, sd=sigma_err))%>%
      mutate(Y=Y-0.2*S)
    
  }else{
    
    # extra phase I data
    lsimXY=simXY%>%
      left_join(data.frame(id=rep(simXY$id, each=6), time=rep(1:6, size)), by="id")
    
    lsimXY=lsimXY%>%
      mutate(Y=sapply(predict(out.model, newdata = lsimXY), rnorm, n=1, sd=sigma_err))
    
  }
  
  
  # long to wide (overwrite)
  simXY=lsimXY%>%
    mutate(time=paste("Y",time,sep=""))%>%
    spread(key = time, value = Y)
  
  
  return(list(simXY))
  
  
}
# result=cal.simXY()
# list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
