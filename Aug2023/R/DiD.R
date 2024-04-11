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
              out.model="Y~X1+X2+X3+X4+X2*X4+X5+S*time+A*time",
              ps.model="S~X1+X2+X3+X4+X3*X4+X5",
              piA.model="A~X1+X2+X3+X4+X5"){
  
   
   
  
  if(method=="DID-EC-OR"){
     
    if(setting %in% c(1,2,3,4)){
      
      ldf=df%>%gather(key = "time", value="Y",Y1:Y4)%>%
        mutate(time=as.numeric(substr(time,2,2)))%>%
        arrange(id, time)%>%
        filter(!(S==1 & A==0 & time>=3)) # exclude RCT control second peiord data
      
      out.model.fit=lm(out.model, data=ldf)
       
      # predicted value for trial subjects
      tmp=ldf%>%filter(S==1)%>%
        mutate(muS1A1=predict(out.model.fit,newdata =mutate(filter(ldf,S==1),S=1,A=1)),
               muS1A0=predict(out.model.fit,newdata =mutate(filter(ldf,S==1),S=1,A=0)),
               muS0A0=predict(out.model.fit,newdata =mutate(filter(ldf,S==1),S=0,A=0)))
      
      delta1=tmp%>%
        filter(time %in% c(1,2))%>%
        mutate(delta=muS1A0-muS0A0)%>%
        group_by(id)%>%
        summarise(delta=mean(delta))%>%
        ungroup()%>%
        summarise(delta=mean(delta))%>%pull()
      
      tau=tmp%>%
        filter(time %in% c(3,4))%>%
        mutate(delta=muS1A1-muS0A0)%>%
        group_by(time)%>%
        summarise(delta=mean(delta)-delta1)%>%pull(delta)
      
    }else{
      ldf=df%>%gather(key = "time", value="Y",Y1:Y6)%>%
        mutate(time=as.numeric(substr(time,2,2)))%>%
        arrange(id, time)%>%
        filter(!(S==1 & A==0 & time>=5)) # exclude RCT control second peiord data
      
      out.model.fit=lm(out.model, data=ldf)
      
      # predicted value for trial subjects
      tmp=ldf%>%filter(S==1)%>%
        mutate(muS1A1=predict(out.model.fit,newdata =mutate(filter(ldf,S==1),S=1,A=1)),
               muS1A0=predict(out.model.fit,newdata =mutate(filter(ldf,S==1),S=1,A=0)),
               muS0A0=predict(out.model.fit,newdata =mutate(filter(ldf,S==1),S=0,A=0)))
      
      delta1=tmp%>%
        filter(time %in% c(1,2,3,4))%>%
        mutate(delta=muS1A0-muS0A0)%>%
        group_by(id)%>%
        summarise(delta=mean(delta))%>%
        ungroup()%>%
        summarise(delta=mean(delta))%>%pull()
      
      tau=tmp%>%
        filter(time %in% c(5,6))%>%
        mutate(delta=muS1A1-muS0A0)%>%
        group_by(time)%>%
        summarise(delta=mean(delta)-delta1)%>%pull(delta)
      
    }
     
      
     
    
  }else if(method=="DID-EC-IPW"){
    
    if(setting %in% c(1,2,3,4)){
      piS=glm(ps.model, data =df, family = "binomial")
      #piA=df%>%filter(S==1)%>%summarise(tmp=mean(A))%>%pull()
      piA=glm(piA.model, data =df[df$S==1,], family = "binomial")
      
      
      temp=df%>%
        mutate(`piAX`=predict(piA,newdata =df,type="response"),
               `piSX`=predict(piS,newdata =df,type="response"),
               rx=`piSX`/(1-`piSX`))%>%
        mutate(w11=1/`piAX`,w10=1/(1-`piAX`),w00=rx)
      
    
      ### create outcomes: obs * T
      Ys=as.matrix(temp[c("Y1","Y2","Y3","Y4")])
      
      potential= (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11)+temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10)+(1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00))*Ys
      mu1=colSums(potential[temp$S==1&temp$A==1,3:4])
      mu0=colSums(potential[temp$S==0,3:4])
      delta1=sum(rowMeans(potential[temp$S==1&temp$A==0,1:2]))-sum(rowMeans(potential[temp$S==0,1:2]))
      
      tau=(mu1-mu0)-delta1
      names(tau)=NULL
    }else{
      piS=glm(ps.model, data =df, family = "binomial")
      #piA=df%>%filter(S==1)%>%summarise(tmp=mean(A))%>%pull()
      piA=glm(piA.model, data =df[df$S==1,], family = "binomial")
      
      temp=df%>%
        mutate(`piAX`=predict(piA,newdata =df,type="response"),
               `piSX`=predict(piS,newdata =df,type="response"),
               rx=`piSX`/(1-`piSX`))%>%
        mutate(w11=1/`piAX`,w10=1/(1-`piAX`),w00=rx)
      
      ### create outcomes: obs * T
      Ys=as.matrix(temp[c("Y1","Y2","Y3","Y4","Y5","Y6")])
      
      potential= (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11)+temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10)+(1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00))*Ys
      mu1=colSums(potential[temp$S==1&temp$A==1,5:6])
      mu0=colSums(potential[temp$S==0,5:6])
      delta1=sum(rowMeans(potential[temp$S==1&temp$A==0,1:4]))-sum(rowMeans(potential[temp$S==0,1:4]))
      
      tau=(mu1-mu0)-delta1
      names(tau)=NULL
    }
     
  }else if(method=="DID-EC-AIPW"){
    
    if(setting %in% c(1,2,3,4)){
      
      ldf=df%>%gather(key = "time", value="Y",Y1:Y4)%>%
        mutate(time=as.numeric(substr(time,2,2)))%>%
        arrange(id, time)%>%
        filter(!(S==1 & A==0 & time>=3)) # exclude RCT control second peiord data
      
      out.model.fit=lm(out.model, data=ldf)

      piS=glm(ps.model, data =df, family = "binomial")
      #piA=df%>%filter(S==1)%>%summarise(tmp=mean(A))%>%pull()
      piA=glm(piA.model, data =df[df$S==1,], family = "binomial")
      
      
      temp=df%>%
        mutate(`piAX`=predict(piA,newdata =df,type="response"),
               `piSX`=predict(piS,newdata =df,type="response"),
               rx=`piSX`/(1-`piSX`))%>%
        mutate(w11=1/`piAX`,w10=1/(1-`piAX`),w00=rx)%>%
        mutate(muS0A01=predict(out.model.fit, newdata =mutate(df,S=0,A=0, time=1)),
               muS0A02=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=2)),
               muS0A03=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=3)),
               muS0A04=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=4)))%>%
        mutate(Y1=Y1-muS0A01,
               Y2=Y2-muS0A02,
               Y3=Y3-muS0A03,
               Y4=Y4-muS0A04)
      
      ### create outcomes: obs * T
      Ys=as.matrix(temp[c("Y1","Y2","Y3","Y4")])
      
      potential= (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11)+temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10)+(1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00))*Ys
      mu1=colSums(potential[temp$S==1&temp$A==1,3:4])
      mu0=colSums(potential[temp$S==0,3:4])
      delta1=sum(rowMeans(potential[temp$S==1&temp$A==0,1:2]))-sum(rowMeans(potential[temp$S==0,1:2]))
      
      tau=(mu1-mu0)-delta1
      names(tau)=NULL
    }else{
      ldf=df%>%gather(key = "time", value="Y",Y1:Y6)%>%
        mutate(time=as.numeric(substr(time,2,2)))%>%
        arrange(id, time)%>%
        filter(!(S==1 & A==0 & time>=5)) # exclude RCT control second peiord data
      
      out.model.fit=lm(out.model, data=ldf)
      
      piS=glm(ps.model, data =df, family = "binomial")
      #piA=df%>%filter(S==1)%>%summarise(tmp=mean(A))%>%pull()
      piA=glm(piA.model, data =df[df$S==1,], family = "binomial")
      
      temp=df%>%
        mutate(`piAX`=predict(piA,newdata =df,type="response"),
               `piSX`=predict(piS,newdata =df,type="response"),
               rx=`piSX`/(1-`piSX`))%>%
        mutate(w11=1/`piAX`,w10=1/(1-`piAX`),w00=rx)%>%
        mutate(muS0A01=predict(out.model.fit, newdata =mutate(df,S=0,A=0, time=1)),
               muS0A02=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=2)),
               muS0A03=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=3)),
               muS0A04=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=4)),
               muS0A05=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=5)),
               muS0A06=predict(out.model.fit,newdata =mutate(df,S=0,A=0, time=6)))%>%
        mutate(Y1=Y1-muS0A01,
               Y2=Y2-muS0A02,
               Y3=Y3-muS0A03,
               Y4=Y4-muS0A04,
               Y5=Y5-muS0A05,
               Y6=Y6-muS0A06)
      
      ### create outcomes: obs * T
      Ys=as.matrix(temp[c("Y1","Y2","Y3","Y4","Y5","Y6")])
      
      potential= (temp$S*temp$A*temp$w11/sum(temp$S*temp$A*temp$w11)+temp$S*(1-temp$A)*temp$w10/sum(temp$S*(1-temp$A)*temp$w10)+(1-temp$S)*temp$w00/sum((1-temp$S)*temp$w00))*Ys
      mu1=colSums(potential[temp$S==1&temp$A==1,5:6])
      mu0=colSums(potential[temp$S==0,5:6])
      delta1=sum(rowMeans(potential[temp$S==1&temp$A==0,1:4]))-sum(rowMeans(potential[temp$S==0,1:4]))
      
      tau=(mu1-mu0)-delta1
      names(tau)=NULL
    }
    
    
  }
   
  
  ####### Use Bootstrap for standard error and confidence intervals
  
    boot.out <- boot(data=df, DiDboot,
                     setting=setting, 
                     ps.model=ps.model,
                     piA.model=piA.model,
                     out.model=out.model,
                     method=method,
                     # parallel = "multicore",
                     # ncpus=4,
                     strata = df$S,
                     R=1000)
  
   
    
  return(list(c(tau[1], sd(boot.out$t[,1]), boot.ci(boot.out,index=1,type=c("perc"))$percent[4:5]),
              c(tau[2], sd(boot.out$t[,2]), boot.ci(boot.out,index=2,type=c("perc"))$percent[4:5])))
            
  
}
 
