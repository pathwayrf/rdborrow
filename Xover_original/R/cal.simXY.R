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
#' result = cal.simXY(setting=2)
#' list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)



cal.simXY = function(path = mydir,
                   setting = 1,
                   size = 300,
                   form_x_s = "SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0",
                   form_x = "SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0+A",
                   form_x_mis = "SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Age_Enrollment^2+Y0+A",
                   form_x_s_mis = "SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Age_Enrollment^2+Y0"){
  
  
  # load the pre-proceed data
  load("Study/output/bene.RData") # patient info
  load("Study/output/outcome.RData") # outcome measures: each row is patient-measure, so one patient has multiple lines corresponding to different measures at different times
  
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
  fit=lm(SMA_Type~1, data=mydat) # 0.3032 SMA Tyoe III
  SMA_Type=sample (c(0,1), size=size, replace=T, prob=c(1-fit$coefficients,fit$coefficients))
  simXY=data.frame(id=1:size,SMA_Type)
  
  mylogit <- glm( Scoliosis ~SMA_Type, data = mydat, family = "binomial")
  pred <- predict(mylogit, newdata = simXY,type='response'  )
  simXY$Scoliosis=sapply(pred, function(prob) sample(c(0,1), size=1,prob=c(1-prob,prob)))
  
  library(nnet) # needed for multinom()
  mylogit <- multinom(  SMN2_Copy_Number~SMA_Type+Scoliosis, data = mydat)
  pred <- predict(mylogit, newdata = simXY, type="probs")
  simXY$SMN2_Copy_Number=sapply(split(pred,row(pred)), function(prob)sample(x=c(2,3,4), size=1, prob=prob))
  
  fit=lm(Age_Enrollment~SMN2_Copy_Number+SMA_Type+Scoliosis, data=mydat)
  mu <- predict(fit, newdata = simXY)
  simXY$Age_Enrollment=sapply(mu,function(x)rnorm(n=1,mean=x,sd=summary(fit)$sigma))
  
  fit=lm(Y0~Age_Enrollment+SMN2_Copy_Number+SMA_Type+Scoliosis, data=mydat)
  mu <- predict(fit, newdata = simXY)
  simXY$Y0=sapply(mu, function(x) rnorm(n=1,mean=x,sd=summary(fit)$sigma))
  
  
  if(setting==1){
    
    m0 <- glm(as.formula(paste("S", "~", form_x_s)), data = mydat, family = "binomial") # propensity of trial participation model
    # outcome models, time specific
    model.list=lapply(1:4, function(x){
      assign(paste("m",x,''),lm(as.formula(paste(paste("Y",x,sep=''), "~", form_x)), data=mydat))
    })
    list2env(setNames(model.list,c("m1","m2","m3","m4")), envir = .GlobalEnv)
    
  }else if(setting==2){
    m0 <- glm(as.formula(paste("S", "~", form_x_s)), data = mydat, family = "binomial")
    # unmeasured confounding with time-invariant effect, all models with observed covariates are correctly specified
    model.list=lapply(1:4, function(x){
      assign(paste("m",x,''),lm(as.formula(paste(paste("Y",x,sep=''), "~", form_x)), data=mydat))
    })
    list2env(setNames(model.list,c("m1","m2","m3","m4")), envir = .GlobalEnv)
    
    #make time-invariant unmeasured confounder : SMA_Type
    models <- list(m2 = m2, m3 = m3, m4 = m4)
    updated_models <- lapply(models, function(model) {
      model$coefficients[2] <- m1$coefficients[2]
      return(model)
    })
    list2env(updated_models, envir = .GlobalEnv)
    
    
    
    
  }else if(setting==3){
    m0 <- glm(as.formula(paste("S", "~", form_x_s_mis)), data = mydat, family = "binomial")
    # unmeasured confounding with time-invariant effect, outcome models with observed covariates are correctly specified
    model.list=lapply(1:4, function(x){
      assign(paste("m",x,''),lm(as.formula(paste(paste("Y",x,sep=''), "~", form_x)), data=mydat))
    })
    list2env(setNames(model.list,c("m1","m2","m3","m4")), envir = .GlobalEnv)
    
    #make time-invariatn unmeasured confounder : SMA_Type
    models <- list(m2 = m2, m3 = m3, m4 = m4)
    updated_models <- lapply(models, function(model) {
      model$coefficients[2] <- m1$coefficients[2]
      return(model)
    })
    list2env(updated_models, envir = .GlobalEnv)
    
    
    
  }else if(setting==4){
    m0 <- glm(as.formula(paste("S", "~", form_x_s)), data = mydat, family = "binomial")
    # unmeasured confounding with time-invariant effect, PS models with observed covariates are correctly specified
    model.list=lapply(1:4, function(x){
      assign(paste("m",x,''),lm(as.formula(paste(paste("Y",x,sep=''), "~", form_x_mis)), data=mydat))
    })
    list2env(setNames(model.list,c("m1","m2","m3","m4")), envir = .GlobalEnv)
    
    #make time-invariatn unmeasured confounder : SMA_Type
    models <- list(m2 = m2, m3 = m3, m4 = m4)
    updated_models <- lapply(models, function(model) {
      model$coefficients[2] <- m1$coefficients[2]
      return(model)
    })
    list2env(updated_models, envir = .GlobalEnv)
    
    
  }else if(setting==5){
    m0 <- glm(as.formula(paste("S", "~", form_x_s_mis)), data = mydat, family = "binomial")
    # unmeasured confounding with time-invariant effect, both models with observed covariates are slightly misspecified
    model.list=lapply(1:4, function(x){
      assign(paste("m",x,''),lm(as.formula(paste(paste("Y",x,sep=''), "~", form_x_mis)), data=mydat))
    })
    list2env(setNames(model.list,c("m1","m2","m3","m4")), envir = .GlobalEnv)
    
    #make time-invariatn unmeasured confounder : SMA_Type
    models <- list(m2 = m2, m3 = m3, m4 = m4)
    updated_models <- lapply(models, function(model) {
      model$coefficients[2] <- m1$coefficients[2]
      return(model)
    })
    list2env(updated_models, envir = .GlobalEnv)
    
  }else if(setting==6){
    m0 <- glm(as.formula(paste("S", "~", form_x_s_mis)), data = mydat, family = "binomial")
    # unmeasured confounding with time-variant effect, both models with observed covariates are slightly misspecified
    # Y is generated with AR(lag=1) potentially depends on previous time point outcome
    
    # time-varying unmeasured confounder: SMA_Type and SMN2_Copy_Number is more promising than Scoliosis, use SMM2_Copy number instead
    m0 <- glm(as.formula(paste("S", "~", form_x_s)), data = mydat, family = "binomial")
    m1=lm(as.formula(paste("Y1", "~", form_x)), data=mydat)
    m2=lm(as.formula(paste("Y2", "~", paste(form_x,"+Y1",sep=""))), data=mydat)  
    m3=lm(as.formula(paste("Y3", "~", paste(form_x,"+Y2",sep=""))), data=mydat)
    m4=lm(as.formula(paste("Y4", "~", paste(form_x,"+Y3",sep=""))), data=mydat)
    
  }
  
  
  if(setting==6){
    m1$coefficients[7]=1
    m2$coefficients[7]=0.5
    m3$coefficients[7]=0.5
    m4$coefficients[7]=0.5
    
    # transmitted via previous outcomes
    m2$coefficients[8]=1
    m3$coefficients[8]=1
    m4$coefficients[8]=1
    # m1$coefficients[7]=1
    # m1$coefficients[7]=1
    # models <- list(m2 = m2, m3 = m3, m4 = m4)
    # updated_models <- lapply(models, function(model) {
    #   model$coefficients[7] <-0.5
    #   # AR coefficient: transmitted via previous outcomes
    #   model$coefficients[8]=1
    #   return(model)
    # })
    # list2env(updated_models, envir = .GlobalEnv)
    
  }else{
    m1$coefficients[7]=1
    m2$coefficients[7]=1.5
    m3$coefficients[7]=2
    m4$coefficients[7]=2.5
  }
   
  
  # Set residual variation  
  sigma_err=1   #1.936
  sigma_b=0.1
  
  # Simulation trial particiaption and treatment assignemnts
  
  # simulate S
  pred <-predict(m0, newdata =simXY,type="response")
  S=sapply(pred, function(prob)sample(c(0,1), size=1, prob=c(1-prob,prob)))
  # simulate A: among S==1 randomization!
  A=S*sample(c(0,1), size,prob=c(1-mean(filter(mydat,S==1)$A),mean(filter(mydat,S==1)$A)), replace = TRUE)
  study=ifelse(S==1 & A==0,"RCT Control Arm",ifelse(S==1 & A==1,"RCT Trt Arm","External Control"))
  simXY=cbind(simXY,study,S,A)
  
  # Simulate outcome Ys
  if(setting==6){
    b=rnorm(size,0,sd=sqrt(sigma_b)) # random effects shared across time
    simXY$Y1=sapply(predict(m1, newdata = simXY)+b, rnorm, n=1, sd=sigma_err)
    simXY$Y2=sapply(predict(m2, newdata = simXY)+b, rnorm, n=1, sd=sigma_err)
    simXY$Y3=sapply(predict(m3, newdata = simXY)+b, rnorm, n=1, sd=sigma_err)
    simXY$Y4=sapply(predict(m4, newdata = simXY)+b, rnorm, n=1, sd=sigma_err)
   
  }else{
    b=rnorm(size,0,sd=sqrt(sigma_b)) # random effects shared across time
    
    Ys.list=lapply(1:4, function(x){
      sapply(predict(get(paste("m",x,sep="")), newdata = simXY)+b, rnorm, n=1, sd=sigma_err)
    })
    list2env(setNames(Ys.list,c("Y1","Y2","Y3","Y4")), envir = .GlobalEnv)
    
    simXY=cbind(simXY,Y1,Y2,Y3,Y4)
  }
  
   
  
  # Get the true ATE in this simulation
   
  if(setting==6){
    ATE=cumsum(c(m1$coefficients[7],m2$coefficients[7],m3$coefficients[7],m4$coefficients[7]))
  }else{
    ATE=sapply(list(m1,m2,m3,m4), function(x) x$coefficients[7])
  }
  
   
  
  
  return(list(simXY,ATE))
  
  
}
# result=cal.simXY()
# list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)

