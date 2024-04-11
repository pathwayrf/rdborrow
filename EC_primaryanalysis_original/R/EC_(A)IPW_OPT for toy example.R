#' Simplified implementation of EC_IPW_OPT() for toy simulation only
#'
#' This function calculate the estimated ATE by the methods EC-IPW with optimal weight as an option, and doubly robust version EC-AIPW with optimal weight as an option.
#'
#' @param df A data frame contains simulated data, returned by cal.simXY()
#' @param wt An initial assignment of weight for external control
#' @param method:
#' "within trial": simple difference in means only using RCT data
#' "EC-IPW-OPT": EC-IPW, if optimal.weight=T, this is EC-IPW-OPT
#' "EC-AIPW-OPT": EC-AIPW, if optimal.weight=T, this is EC-AIPW-OPT
#' @param optimal.weight Use optimal weight or pre-specified weight?
#' @param form_x Assumed working model for nuisance components, such as propensity score model and outcome models
#' @return A list contains: estimated ATE, SE, weihgt used 
#' @examples
#' mysimXY=toy.simXY()
#' EC_IPW_OPT_toy(df=mysimXY,wt=0,method="EC-IPW-OPT")
#' EC_IPW_OPT_toy(df=mysimXY,wt=0,method="EC-AIPW-OPT")
#' EC_IPW_OPT_toy(df=mysimXY,wt=0,method="within trial")
#' 

EC_IPW_OPT_toy<-function(df,
                     wt=0,
                     method="EC-IPW-OPT",
                     optimal.weight=T,
                     form_x="X"){

  n=sum(df$S) # RCT sample size
  m=sum(1-df$S) # external control sample size
  pi.S=n/(n+m) # marginal propensity of trial particiaption
  

  #============= within trial ===============
  if(method=="within trial"){
    temp=df%>%
      filter(S==1)%>%
      mutate(`piA`=sum(A[S==1])/n)%>%
      mutate(w11=`piA`,w10=1-`piA`)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y")])

    potential= (temp$A/temp$w11+(1-temp$A)/temp$w10)*Ys
    mu1=sum(potential[temp$S==1&temp$A==1])/n
    mu10=sum(potential[temp$S==1&temp$A==0])/n
     

    # variance

    # bread
    A=diag(c(-1,-1))

    # meet
    phi1=temp$A*(Ys-mu1)/temp$`piA`     # influence from rct treated
    phi2=(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)     # influence from rct control
     
    ## big Phi should be n by 4+2+p, observations of same phi in the same column
    Phi=cbind(phi1,phi2)

    B=(t(Phi)%*%Phi)/(n)

    sigma=solve(A)%*%B%*%t(solve(A))
  
    tau=mu1-mu10

    coef.mat=c(1,-1)

    sd.tau=sqrt(coef.mat%*%sigma%*%coef.mat/(n))


  } # =============== EC-AIPW-OPT ==============
  else if(method=="EC-AIPW-OPT"){
    # propensity score model
    piS.model=glm(as.formula(paste("S","~",form_x)) , data =df, family = "binomial")
    # outcome regression model
    Y0.model=lm(as.formula(paste("Y","~",form_x)), data = filter(df,A==0) )
    
    
    
    temp=df%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model,newdata =df,type="response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`),
             `Y0`=predict(Y0.model,newdata =df),
             Y=Y-`Y0`)%>%
      mutate(w11=`piA`,w10=1-`piA`,w00=rx)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y")])

    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1=sum(potential[temp$S==1&temp$A==1,])/n
    mu10=sum(potential[temp$S==1&temp$A==0,])/n
    mu00=sum(potential[temp$S==0,])/m


    # variance

    # bread
    A=diag(c(-1,-1,-mean((1-temp$S)*temp$w00/(1-temp$piS))))

    # meet
    phi1=temp$S*temp$A*(Ys-mu1)/temp$`piA`/temp$piS    # influence from rct treated
    phi2=temp$S*(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)/temp$piS    # influence from rct control
    phi3=(1-temp$S)*(Ys-mu00)*temp$w00/(1-temp$piS)

    ## big Phi should be n by 4+2+p, observations of same phi in the same column
    Phi=cbind(phi1,phi2,phi3 )

    B=(t(Phi)%*%Phi)/(n+m)

    sigma=solve(A)%*%B%*%t(solve(A))

    # optimal weight =
    #w.opt=sigma[4,4]/(sigma[4,4]+sigma[6,6])
    w.opt=(sigma[2,2])/(sigma[2,2]+sigma[3,3]+(n+m)*(mu10-mu00)^2)
    if(optimal.weight==T){wt=w.opt}


    mu0=(1-wt)*mu10+wt*mu00
    tau=mu1-mu0

    coef.mat=c(1,-(1-wt),-wt)
    sd.tau=sqrt(coef.mat%*%sigma%*%coef.mat/(n+m))




  } # ============== EC-IPW-OPT ==================
  else if(method=="EC-IPW-OPT"){

    piS.model=glm(as.formula(paste("S","~",form_x)) , data =df, family = "binomial")

    # estimate ATE
    temp=df%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model,newdata =df,type="response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`))%>%
      mutate(w11=`piA`,w10=1-`piA`,w00=rx)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y")])

    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1=sum(potential[temp$S==1&temp$A==1,])/n
    mu10=sum(potential[temp$S==1&temp$A==0,])/n
    mu00=sum(potential[temp$S==0,])/m


    # variance

    # bread
    A34=t((as.vector((1-temp$S)*temp$`piSX`/(temp$`piS`*(1-temp$`piSX`)))*(Ys-mu00)))%*%model.matrix(piS.model)/(n+m)
    piS.beta=t(model.matrix(piS.model))%*%diag(-temp$`piSX`*(1-temp$`piSX`))%*%model.matrix(piS.model)/(n+m)
    A44=piS.beta

    A<-matrix(0,nrow=5,ncol=5)
    A[1,1]=-1
    A[2,2]=-1
    A[3,3]=-mean((1-temp$S)*temp$w00/(1-temp$piS))
    A[3,4:5]=A34
    A[4:5,4:5]=A44


    # meet
    phi1=temp$S*temp$A*(Ys-mu1)/temp$`piA`/temp$piS    # influence from rct treated
    phi2=temp$S*(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)/temp$piS    # influence from rct control
    phi3=(1-temp$S)*(Ys-mu00)*temp$w00/(1-temp$piS)
    phi.piS=(temp$S-temp$`piSX`)*model.matrix(piS.model)

    ## big Phi should be n by 4+2+p, observations of same phi in the same column
    Phi=cbind(phi1,phi2,phi3,phi.piS )

    B=(t(Phi)%*%Phi)/(n+m)

    sigma=solve(A)%*%B%*%t(solve(A))

    # optimal weight =
    #w.opt=sigma[4,4]/(sigma[4,4]+sigma[6,6])
    w.opt=(sigma[2,2]+sigma[1,3]-sigma[1,2]-sigma[2,3])/(sigma[2,2]+sigma[3,3]-2*sigma[2,3]+(n+m)*(mu10-mu00)^2)
    if(optimal.weight==T){wt=w.opt}


    mu0=(1-wt)*mu10+wt*mu00
    tau=mu1-mu0

    coef.mat=c(1,-(1-wt),-wt,0,0)
    sd.tau=sqrt(coef.mat%*%sigma%*%coef.mat/(n+m))



  }else if(method=="EC-AIPW-OPT"){
    piS.model=glm(as.formula(paste("S","~","X")) , data =df, family = "binomial")
    Y0.model=lm(as.formula(paste("Y","~","X")), data = filter(df,A==0) )

    # estimate ATE
    temp=df%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model,newdata =df,type="response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`),
             `Y0`=predict(Y0.model,newdata =df),
             Y=Y-`Y0`)%>%
      mutate(w11=`piA`,w10=1-`piA`,w00=rx)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y")])


    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1=sum(potential[temp$S==1&temp$A==1,])/n
    mu10=sum(potential[temp$S==1&temp$A==0,])/n
    mu00=sum(potential[temp$S==0,])/m

    # variance

    # bread
    A34=t((as.vector((1-temp$S)*temp$`piSX`/(temp$`piS`*(1-temp$`piSX`)))*(Ys-mu00)))%*%model.matrix(piS.model)/(n+m)
    piS.beta=t(model.matrix(piS.model))%*%diag(-temp$`piSX`*(1-temp$`piSX`))%*%model.matrix(piS.model)/(n+m)
    A44=piS.beta

    phi1.gamma=as.vector(-temp$S*temp$A/(temp$`piS`*(temp$`piA`)))%*%model.matrix(piS.model)/(n+m) # gamma
    phi2.gamma=as.vector(-temp$S*(1-temp$A)/(temp$`piS`*(1-temp$`piA`)))%*%model.matrix(piS.model)/(n+m) # gamma
    phi3.gamma=as.vector(-(1-temp$S)*temp$rx/(1-temp$`piS`))%*%model.matrix(piS.model)/(n+m) # gamma
    Y0.gamma=-t(model.matrix(piS.model))%*%diag((1-temp$A)/(1-mean(temp$A)))%*%model.matrix(piS.model)/(n+m)
    A55=Y0.gamma

    A<-matrix(0,nrow=7,ncol=7)
    A[1,1]=-1
    A[2,2]=-1
    A[3,3]=-mean((1-temp$S)*temp$w00/(1-temp$piS))
    A[3,4:5]=A34
    A[4:5,4:5]=A44
    A[6:7,6:7]=A55


    A[1,6:7]=phi1.gamma
    A[2,6:7]=phi2.gamma
    A[3,6:7]=phi3.gamma

    # meet
    phi1=temp$S*temp$A*(Ys-mu1)/temp$`piA`/temp$piS    # influence from rct treated
    phi2=temp$S*(1-temp$A)*(Ys-mu10)/(1-temp$`piA`)/temp$piS    # influence from rct control
    phi3=(1-temp$S)*(Ys-mu00)*temp$w00/(1-temp$piS)
    phi.piS=(temp$S-temp$`piSX`)*model.matrix(piS.model)
    phi.Y0=((1-temp$A)/(1-mean(temp$A)))*(Ys[,1]*model.matrix(piS.model))

    ## big Phi should be n by 4+2+p, observations of same phi in the same column
    Phi=cbind(phi1,phi2,phi3,phi.piS,phi.Y0 )

    B=(t(Phi)%*%Phi)/(n+m)

    sigma=solve(A)%*%B%*%t(solve(A))

    # optimal weight =
    #w.opt=sigma[4,4]/(sigma[4,4]+sigma[6,6])
    w.opt=(sigma[2,2]+sigma[1,3]-sigma[1,2]-sigma[2,3])/(sigma[2,2]+sigma[3,3]-2*sigma[2,3]+(n+m)*(mu10-mu00)^2)
    if(optimal.weight==T){wt=w.opt}


    mu0=(1-wt)*mu10+wt*mu00
    tau=mu1-mu0
    coef.mat=c(1,-(1-wt), -wt,0,0,0,0)

    sd.tau=sqrt(coef.mat%*%sigma%*%coef.mat/(n+m))



  }

    return(c(tau, sd.tau, wt))



}

