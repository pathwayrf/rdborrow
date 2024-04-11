#' Implement the EC-IPW(-OPT) and EC-AIPW(-OPT) in the manuscript "Causal Estimators for Incorporating External Controls in Randomized Trials with Longitudinal Outcomes"
#'
#' This function calculate the estimated ATE by the methods EC-IPW with optimal weight as an option, and doubly robust version EC-AIPW with optimal weight as an option.
#'
#' @param df A data frame contains simulated data, returned by cal.simXY()
#' @param wt An initial assignment of weight for external control
#' @param method:
#' "within trial": simple difference in means only using RCT data
#' "EC-IPW-OPT": EC-IPW, if optimal.weight=T, this is EC-IPW-OPT
#' "EC-AIPW-OPT": EC-AIPW, if optimal.weight=T, this is EC-AIPW-OPT
#' @param setting The setting number
#' Setting 1: all confounder correctly specified
#' Setting 2: outcome model misspecifed
#' Setting 3: propensity model misspecifed
#' Setting 4: both propensity model and outcome models misspecified
#' Setting 5: Exists unmeasured confounder: SMM2_Copy number as pseudo omitted confounder
#' @param Bootstrap Also bootstrap SE?
#' @param optimal.weight Use optimal weight or pre-specified weight?
#' @param form_x Assumed working model for nuisance components, such as propensity score model and outcome models
#' @param form_x_omit Only used when setting==5 where one covariates is intentionally treated as unmeasured
#' @return A list contains: estimated ATE, SE, weihgt used, SE by Bootstrap and a 95% confidence interval for primary endpoint (only when Bootstrap=TRUE)
#' @examples
#' result=cal.simXY()
#' list2env(setNames(result,c("df","ATE")), envir = .GlobalEnv)
#' EC_IPW_OPT(df=mysimXY,wt=0,method="EC-IPW-OPT",setting=1,Bootstrap=F,optimal.weight=T)
#' EC_IPW_OPT(df=mysimXY,wt=0,method="EC-AIPW-OPT",setting=1,Bootstrap=F,optimal.weight=T)
#' EC_IPW_OPT(df=mysimXY,wt=0,method="EC-IPW-OPT",setting=1,Bootstrap=T,optimal.weight=T)
#' EC_IPW_OPT(df=mysimXY,wt=0,method="EC-AIPW-OPT",setting=1,Bootstrap=T,optimal.weight=T)

EC_IPW_OPT_old <- function(df,
                     wt=0,
                     method="EC-IPW-OPT",
                     setting=1,
                     Bootstrap=F,
                     optimal.weight=T,
                     form_x="SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0",
                     form_x_omit="SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0"){


  n=sum(df$S) # RCT sample size
  m=sum(1-df$S) # external control sample size
  pi.S=n/(n+m) # marginal propensity of trial particiaption

  # model form for propensity score or outcome regression model (linear in baseline covariates), omit SMN2_Copy_Number in setting 5
  if(setting==5){form_x=form_x_omit} 
  
  if(method=="within trial"){
    # estimate ATE
    temp=df%>%
      filter(S==1)%>%
      mutate(`piA`=sum(A)/n)%>%
      mutate(w11=`piA`,w10=1-`piA`)

    ### create outcomes: obs by T
    # Ys=as.matrix(temp[c("Y1","Y2")])
    Ys=as.matrix(temp[c("Y1","Y2")])

    potential= (temp$A/temp$w11+(1-temp$A)/temp$w10)*Ys
    mu1=colSums(potential[temp$A==1,])/n
    mu0=colSums(potential[temp$A==0,])/n

    tau=mu1-mu0

    # asymptotic variance of ATE
    ## bread
    A=diag(rep(-1,2*dim(Ys)[2])  )
    ## meet
    phi1=temp$A*(Ys-mu1)/temp$`piA`  # influence from rct treated
    phi2=(1-temp$A)*(Ys-mu0)/(1-temp$`piA`)  # influence from rct control
    ## big Phi should be n by T
    Phi=cbind(phi1,phi2 )
    B=(t(Phi)%*%Phi)/n
    ## sandwich
    sigma=solve(A)%*%B%*%t(solve(A))

    # ATE as linear comb of parameters
    coef.mat=rbind(c(1,0, -1, 0),c(0,1, 0,-1))
    # get standard error using the linear comb of var-cov matrix
    sd.tau=sqrt(diag(coef.mat%*%sigma%*%t(coef.mat)/(n)))

  }else if(method=="EC-IPW-OPT"){
    # propensity score model
    piS.model=glm(as.formula(paste("S","~",form_x)) , data =df, family = "binomial")
    # estimate ATE
    temp=df%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model,newdata =df,type="response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`))%>%
      mutate(w11=`piA`,w10=1-`piA`,w00=rx)

    ### create outcomes: obs * T
    # Ys=as.matrix(temp[c("Y1","Y2")])
    Ys=as.matrix(temp[c("Y1","Y2")])

    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1=colSums(potential[temp$S==1&temp$A==1,])/sum(temp$S*temp$A/temp$w11)
    mu10=colSums(potential[temp$S==1&temp$A==0,])/sum(temp$S*(1-temp$A)/temp$w10)
    mu00=colSums(potential[temp$S==0,])/(sum((1-temp$S)*temp$w00))


    # variance

    ## bread
    A11=diag(rep(-1,dim(Ys)[2]))
    A22=diag(rep(-1,dim(Ys)[2]))
    A33=diag(rep(-mean((1-temp$S)*temp$w00/(1-temp$piS)),dim(Ys)[2]))
    A34=t((as.vector((1-temp$S)*temp$`piSX`/(temp$`piS`*(1-temp$`piSX`)))*(Ys-mu00)))%*%model.matrix(piS.model)/(n+m)
    piS.beta=t(model.matrix(piS.model))%*%diag(-temp$`piSX`*(1-temp$`piSX`))%*%model.matrix(piS.model)/(n+m)
    A44=piS.beta

    n1=dim(A11)[1];n2=dim(A22)[1];n3=dim(A33)[1];n4=dim(A44)[1]

    A<-matrix(0,nrow=n1+n2+n3+n4,ncol=n1+n2+n3+n4)
    A[1:n1,1:n1]=A11
    A[(n1+1):(n1+n2),(n1+1):(n1+n2)]=A22
    A[(n1+n2+1):(n1+n2+n3),(n1+n2+1):(n1+n2+n3)]=A33
    A[(n1+n2+n3+1):(n1+n2+n3+n4),(n1+n2+n3+1):(n1+n2+n3+n4)]=A44

    A[(n1+n2+1):(n1+n2+n3),(n1+n2+n3+1):(n1+n2+n3+n4)]=A34


    ## meet
    phi1=temp$S*temp$A*(Ys-mu1)/temp$`piA`/temp$piS    # influence from rct treated
    phi2=temp$S*(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)/temp$piS    # influence from rct control
    phi3=(1-temp$S)*(Ys-mu00)*temp$w00/(1-temp$piS)
    phi.piS=(temp$S-temp$`piSX`)*model.matrix(piS.model)

    ## big Phi should be n by T+T+p,
    Phi=cbind(phi1,phi2,phi3,phi.piS )

    B=(t(Phi)%*%Phi)/(n+m)

    ## sandwich
    sigma=solve(A)%*%B%*%t(solve(A))

    ## Optimal weight as proposed in manuscript
    fit1=lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==1&A==0) )
    sigma10=(summary(fit1)$sigma)**2
    num=sum(temp$S*(1-temp$A)/temp$w10^2/(sum(temp$S*(1-temp$A)/temp$w10))^2)
    #+ sum(temp$S*(1-temp$A))*var(temp$S*(1-temp$A)*predict(fit1,newdata =temp)/temp$w10/sum(temp$S*(1-temp$A)/temp$w10))

    fit0=lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==0) )
    sigma00=(summary(fit0)$sigma)**2
    denom=sum((1-temp$S)*temp$w00^2/(sum((1-temp$S)*temp$w00))^2)
    #+ sum((1-temp$S))*var((1-temp$S)*predict(fit0,newdata =temp)*temp$w00/sum((1-temp$S)*temp$w00))
    w.opt=num/(num+denom)
    # only use the optimal weight if we want
    if(optimal.weight==T){wt=w.opt}

    # final hybrid estimate as combination of RCT control and external control
    mu0=(1-wt)*mu10+wt*mu00
    tau=mu1-mu0
    # ATE as linear comb of parameters
    coef.mat=matrix(c(1,0, -(1-wt), 0, -wt,0, rep(0,dim(piS.beta)[2]), 0,1, 0,-(1-wt), 0, -wt, rep(0,dim(piS.beta)[2])),byrow = T,nrow = 2)
    # get standard error using the linear comb of var-cov matrix
    sd.tau=sqrt(diag(coef.mat%*%sigma%*%t(coef.mat)/(n+m)))


  }else if(method=="EC-AIPW-OPT"){
    # propensity score model
    piS.model=glm(as.formula(paste("S","~",form_x)), data =df, family = "binomial")
    # outcome regression model
    Y0.model=list(lm(as.formula(paste("Y1","~",form_x)), data = filter(df,A==0) ),
                  lm(as.formula(paste("Y2","~",form_x)), data = filter(df,A==0) ))


    # estimate ATE
    temp=df%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model, newdata = df, type = "response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`),
             `Y1_0` = predict(Y0.model[[1]], newdata = df),
             `Y2_0` = predict(Y0.model[[2]], newdata = df),
             Y1 = Y1 - `Y1_0`,
             Y2 = Y2 - `Y2_0`) %>%
      mutate(w11=`piA`, w10=1-`piA`, w00=rx)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y1","Y2")])


    potential = (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1 = colSums(potential[temp$S==1&temp$A==1,])/n
    mu10 = colSums(potential[temp$S==1&temp$A==0,])/n
    mu00 = colSums(potential[temp$S==0,])/(sum((1-temp$S)*temp$w00))



    # variance

    ## bread
    A11 = diag(rep(-1,dim(Ys)[2]))
    A22 = diag(rep(-1,dim(Ys)[2]))
    A33 = diag(rep(-mean((1-temp$S)*temp$w00/(1-temp$piS)),dim(Ys)[2]))
    A34 = t((as.vector((1-temp$S)*temp$`piSX`/(temp$`piS`*(1-temp$`piSX`)))*(Ys-mu00)))%*%model.matrix(piS.model)/(n+m)
    piS.beta = t(model.matrix(piS.model))%*%diag(-temp$`piSX`*(1-temp$`piSX`))%*%model.matrix(piS.model)/(n+m)
    A44 = piS.beta

    phi1.gamma=as.vector(-temp$S*temp$A/(temp$`piS`*(temp$`piA`)))%*%model.matrix(piS.model)/(n+m) # gamma
    phi2.gamma=as.vector(-temp$S*(1-temp$A)/(temp$`piS`*(1-temp$`piA`)))%*%model.matrix(piS.model)/(n+m) # gamma
    phi3.gamma=as.vector(-(1-temp$S)*temp$rx/(1-temp$`piS`))%*%model.matrix(piS.model)/(n+m) # gamma
    Y0.gamma=-t(model.matrix(piS.model))%*%diag((1-temp$A)/(1-mean(temp$A)))%*%model.matrix(piS.model)/(n+m)
    A55=Y0.gamma

    n1=dim(A11)[1];n2=dim(A22)[1];n3=dim(A33)[1];n4=dim(A44)[1];n5=2*dim(A55)[1];

    A<-matrix(0,nrow=n1+n2+n3+n4+n5,ncol=n1+n2+n3+n4+n5)
    A[1:n1,1:n1]=A11
    A[(n1+1):(n1+n2),(n1+1):(n1+n2)]=A22
    A[(n1+n2+1):(n1+n2+n3),(n1+n2+1):(n1+n2+n3)]=A33
    A[(n1+n2+n3+1):(n1+n2+n3+n4),(n1+n2+n3+1):(n1+n2+n3+n4)]=A44
    A[(n1+n2+n3+n4+1):(n1+n2+n3+n4+dim(A55)[1]),(n1+n2+n3+n4+1):(n1+n2+n3+n4+dim(A55)[1])]=A55
    A[(n1+n2+n3+n4+dim(A55)[1]+1):(n1+n2+n3+n4+n5),(n1+n2+n3+n4+dim(A55)[1]+1):(n1+n2+n3+n4+n5)]=A55

    A[(1):(n1),(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]=matrix(c(phi1.gamma,rep(0,length(phi1.gamma)),rep(0,length(phi1.gamma)),phi1.gamma),nrow=2,byrow = T)
    A[(n1+1):(n1+n2),(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]=matrix(c(phi2.gamma,rep(0,length(phi2.gamma)),rep(0,length(phi2.gamma)),phi2.gamma),nrow=2,byrow = T)
    A[(n1+n2+1):(n1+n2+n3),(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]=matrix(c(phi3.gamma,rep(0,length(phi3.gamma)),rep(0,length(phi3.gamma)),phi3.gamma),nrow=2,byrow = T)
    A[(n1+n2+1):(n1+n2+n3),(n1+n2+n3+1):(n1+n2+n3+n4)]=A34

    ## meet
    phi1=temp$S*temp$A*(Ys-mu1)/temp$`piA`/temp$piS    # influence from rct treated
    phi2=temp$S*(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)/temp$piS    # influence from rct control
    phi3=(1-temp$S)*(Ys-mu00)*temp$w00/(1-temp$piS)
    phi.piS=(temp$S-temp$`piSX`)*model.matrix(piS.model)
    phi.Y0=cbind(((1-temp$A)/(1-mean(temp$A)))*(Ys[,1]*model.matrix(piS.model)),
                 ((1-temp$A)/(1-mean(temp$A)))*(Ys[,2]*model.matrix(piS.model)))

    ## big Phi should be n by 4+2+p, observations of same phi in the same column
    Phi=cbind(phi1,phi2,phi3,phi.piS,phi.Y0 )

    B=(t(Phi)%*%Phi)/(n+m)
    ## sandwich
    sigma=solve(A)%*%B%*%t(solve(A))

    ## Optimal weight as proposed in manuscript
    fit1=lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==1&A==0) )
    sigma10=(summary(fit1)$sigma)**2
    num=sum(temp$S*(1-temp$A)/temp$w10^2/(sum(temp$S*(1-temp$A)/temp$w10))^2)
    #+ sum(temp$S*(1-temp$A))*var(temp$S*(1-temp$A)*predict(fit1,newdata =temp)/temp$w10/sum(temp$S*(1-temp$A)/temp$w10))

    fit0 = lm(as.formula(paste("Y2","~",form_x)), data = filter(temp,S==0) )
    sigma00 = (summary(fit0)$sigma)**2
    denom = sum((1-temp$S)*temp$w00^2/(sum((1-temp$S)*temp$w00))^2)
    #+ sum((1-temp$S))*var((1-temp$S)*predict(fit0,newdata =temp)*temp$w00/sum((1-temp$S)*temp$w00))
    w.opt=num/(num+denom)
    if(optimal.weight==T){wt=w.opt}

    # final hybrid estimate as combination of RCT control and external control
    mu0=(1-wt)*mu10+wt*mu00
    tau=mu1-mu0

    # ATE as linear comb of parameters
    coef.mat=rbind(c(1,0, -(1-wt), 0, -wt,0, rep(0,dim(piS.beta)[2]+2*dim(Y0.gamma)[2])),
                   c(0,1, 0,-(1-wt), 0, -wt, rep(0,dim(piS.beta)[2]+2*dim(Y0.gamma)[2])))
    # get standard error using the linear comb of var-cov matrix
    sd.tau=sqrt(diag(coef.mat%*%sigma%*%t(coef.mat)/(n+m)))


  }

  # If use Bootstrap to construct confidence interval, include bootstrap standard error using bias-corrected and accelerated method, add it as last elememnt in the returned list
  if(Bootstrap==T){
    library(boot)
    boot.out <- boot(data=df, statistic=EC_IPW_OPT_boot, wt=wt, method=method, setting=setting,
                     R=3000)
    ci=boot.ci(boot.out,
               type = c("bca"), index=2)

    sd.boot=sqrt(diag(var(boot.out$t)))
    return(list(tau, sd.tau, wt, sd.boot,ci$bca[4:5] ))
  }else{
    return(list(tau, sd.tau, wt))
  }


}
