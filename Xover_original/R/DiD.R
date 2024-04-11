#' Implement the 3 DID methods: DID-EC-OR, DID-EC-IPW,and DID-EC-AIPW, in the manuscript  
#'
#' This function calculate the estimated ATE by 3 DID methods and Bootstrap CI.
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
#' @param form_x Assumed working model for nuisance components, such as propensity score model and outcome models
#' @param form_x_omit used when setting=2,3,4,5,6 where one covariates is intentionally treated as unmeasured
#' @return A list contains: estimated ATE, SE, weihgt used, SE by Bootstrap and a 95% confidence interval for primary endpoint (only when Bootstrap=TRUE)
#' @examples
#' result=cal.simXY()
#' list2env(setNames(result,c("df","ATE")), envir = .GlobalEnv)
#' DiD()
#' DiD(method="DID-EC-IPW")
#' DiD(method="DID-EC-AIPW")
 
DiD<-function(df=mysimXY,
              setting=1,
              method="DID-EC-OR",
              form_x="SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0",
              form_x_omit="SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0"){
  
#  piS_data = 
  
  if(method=="DID-EC-OR"){
    if(setting==1){
      # external outcome model
      model.list=lapply(1:4, function(x){
        assign(paste("m.ext",x,''),lm(as.formula(paste(paste("y",x,sep=''), "~", form_x)), data=filter(df,S==0)))
      })
      list2env(setNames(model.list,c("m.ext1","m.ext2","m.ext3","m.ext4")), envir = .GlobalEnv)
      
      # rct outcome model
      m.rct1=lm(as.formula(paste("y1","~",form_x)),data = filter(df,S==1&A==0) )
      m.rct2=lm(as.formula(paste("y2","~",form_x)),data = filter(df,S==1&A==0) )
      m.rct3=lm(as.formula(paste("y3","~",form_x)),data = filter(df,S==1&A==1) )
      m.rct4=lm(as.formula(paste("y4","~",form_x)),data = filter(df,S==1&A==1) )
      
    }else{
      # external outcome model
      model.list=lapply(1:4, function(x){
        assign(paste("m.ext",x,''),lm(as.formula(paste(paste("y",x,sep=''), "~", form_x_omit)), data=filter(df,S==0)))
      })
      list2env(setNames(model.list,c("m.ext1","m.ext2","m.ext3","m.ext4")), envir = .GlobalEnv)
      
      # rct outcome model
      m.rct1=lm(as.formula(paste("y1","~",form_x_omit)),data = filter(df,S==1&A==0) )
      m.rct2=lm(as.formula(paste("y2","~",form_x_omit)),data = filter(df,S==1&A==0) )
      m.rct3=lm(as.formula(paste("y3","~",form_x_omit)),data = filter(df,S==1&A==1) )
      m.rct4=lm(as.formula(paste("y4","~",form_x_omit)),data = filter(df,S==1&A==1) )
      
      
    }
    
    # predicted value for trial subjects
    df %>% filter(S==1)%>%
      mutate(muo1=predict(m.ext1,newdata =filter(df,S==1)),
             muo2=predict(m.ext2,newdata =filter(df,S==1)),
             muo3=predict(m.ext3,newdata =filter(df,S==1)),
             muo4=predict(m.ext4,newdata =filter(df,S==1)),
             mue1=predict(m.rct1,newdata =filter(df,S==1)),
             mue2=predict(m.rct2,newdata =filter(df,S==1)),
             mue3=predict(m.rct3,newdata =filter(df,S==1)),
             mue4=predict(m.rct4,newdata =filter(df,S==1)))%>%
      mutate(`tau(t=3)`=mue3-muo3-((mue1-muo1)+(mue2-muo2))/2,
             `tau(t=4)`=mue4-muo4-((mue1-muo1)+(mue2-muo2))/2)%>%
      summarize(`tau(t=3)`=mean(`tau(t=3)`),
                `tau(t=4)`=mean(`tau(t=4)`)
      )%>%
      ungroup()->tau
    
    tau=unlist(tau)
    
  }else if(method == "DID-EC-IPW"){
    
    if(setting==1){
      piS = glm(as.formula(paste("S","~",form_x)), data = df, family = "binomial")
      piA = glm(as.formula(paste("A","~",form_x)), data = filter(df,S==1), family = "binomial")
    }else{
      piS = glm(as.formula(paste("S","~",form_x_omit)), data =df, family = "binomial")
      piA = glm(as.formula(paste("A","~",form_x_omit)), data =filter(df,S==1), family = "binomial")
    }
    
    temp = df %>%
      mutate(`piAX` = predict(piA,newdata = filter(df), type = "response"), #sum(A)/sum(S)
             `piSX` = predict(piS,newdata = filter(df), type = "response"),
#             rx=`piSX`/(1-`piSX`)) %>%
             rx = piSX * (1 - piS)/(1 - piSX)/piS) %>%
      mutate(w11=1/`piAX`, w10=1/(1-`piAX`), w00=rx)
    
    ### create outcomes: obs * T
    Ys = as.matrix(temp[c("y1","y2","y3","y4")])
    
    potential = (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11)+temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10)+(1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00))*Ys
    mu1 = colSums(potential[temp$S==1 & temp$A==1, 3:4])
    mu0 = colSums(potential[temp$S==0, 3:4])
    bias = sum(rowMeans(potential[temp$S==1 & temp$A==0,1:2])) - sum(rowMeans(potential[temp$S==0,1:2]))
    
    print(c(mu1, mu0, bias))
    
    tau = mu1 - mu0 - bias
    names(tau) = c("tau(t=3)","tau(t=4)")
     
  }else if(method=="DID-EC-AIPW"){
    
    if(setting==1){
      piS = glm(as.formula(paste("S","~",form_x)), data =df, family = "binomial")
      piA = glm(as.formula(paste("A","~",form_x)), data =filter(df,S==1), family = "binomial")
      # external outcome model
      model.list=lapply(1:4, function(x){
        assign(paste("m.ext",x,''),lm(as.formula(paste(paste("y",x,sep=''), "~", form_x)), data=filter(df,S==0)))
      })
      list2env(setNames(model.list,c("m.ext1", "m.ext2", "m.ext3", "m.ext4")), envir = .GlobalEnv)
      
    }else{
      piS=glm(as.formula(paste("S","~",form_x_omit)), data =df, family = "binomial")
      piA=glm(as.formula(paste("A","~",form_x_omit)), data =filter(df, S==1), family = "binomial")
      # external outcome model
      model.list=lapply(1:4, function(x){
        assign(paste("m.ext",x,''),lm(as.formula(paste(paste("y",x,sep=''), "~", form_x_omit)), data=filter(df,S==0)))
      })
      list2env(setNames(model.list,c("m.ext1", "m.ext2", "m.ext3", "m.ext4")), envir = .GlobalEnv)
      
    }
     
     
    
    temp=df%>%
      mutate(
        piAX = predict(piA, newdata = filter(df), type="response"), #sum(A)/sum(S)
        piSX = predict(piS, newdata = filter(df), type="response"),
        rx = piSX/(1 - piSX),
        muo1=predict(m.ext1, newdata = filter(df)),
        muo2=predict(m.ext2, newdata = filter(df)),
        muo3=predict(m.ext3, newdata = filter(df)),
        muo4=predict(m.ext4, newdata = filter(df))) %>%
      mutate(w11 = 1/`piAX`, w10 = 1/(1-`piAX`), w00 = rx)%>%
      mutate(y1=y1-muo1,
             y2=y2-muo2,
             y3=y3-muo3,
             y4=y4-muo4)
    
    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("y1","y2","y3","y4")])
    
    potential = (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11)+temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10)+(1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00))*Ys
    mu1 = colSums(potential[temp$S==1&temp$A==1,3:4])
    mu0 = colSums(potential[temp$S==0,3:4])
    bias = sum(rowMeans(potential[temp$S==1&temp$A==0, 1:2]))-sum(rowMeans(potential[temp$S==0, 1:2]))
    
    
    tau = mu1-mu0-bias
    names(tau) = c("tau(t=3)","tau(t=4)")
    
    
  }
   
  
  ####### Use Bootstrap for standard error and confidence intervals
  
    boot.out <- boot(data=df, DiDboot,
                     setting = setting, 
                     method = method,
                     form_x = form_x,
                     form_x_omit = form_x_omit,
                     # parallel = "multicore",
                     # ncpus=4,
                     strata = df$S,
                     R = 1000)
  
  
  return(list(c(tau[1], sd(boot.out$t[,1]), boot.ci(boot.out,index=1,type=c("perc"))$percent[4:5]),
              c(tau[2], sd(boot.out$t[,2]), boot.ci(boot.out,index=2,type=c("perc"))$percent[4:5])))
            
  
}
 
