#' Bootstrap the EC-IPW(-OPT) and EC-AIPW(-OPT) in the manuscript "Causal Estimators for Incorporating External Controls in Randomized Trials with Longitudinal Outcomes"
#'
#' This function calculate the bootstrap SE for EC-IPW(-OPT) and EC-AIPW(-OPT) and should be used inside of EC_IPW_OPT() when Bootstrap=TRUE.
#'
#' @param data A data frame contains simulated data used inside EC_IPW_OPT()
#' @param wt Weight used inside EC_IPW_OPT()
#' @param method:
#' "within trial": simple difference in means only using RCT data
#' "EC-IPW-OPT": EC-IPW, if optimal.weight=T, this is EC-IPW-OPT
#' "EC-AIPW-OPT": EC-AIPW, if optimal.weight=T, this is EC-AIPW-OPT
#' @return A vector ATEs, one for each time point.
#' @examples
#' boot.out <- boot(data=df, statistic=EC_IPW_OPT_boot, wt=wt, method=method, setting=setting, R=3000)
#'
#'
#'
EC_IPW_OPT_boot<-function(data, indices, wt, method, setting){

  # model form for propensity score or outcome regression model (linear in baseline covariates), omit SMN2_Copy_Number in setting 5
  form_x=ifelse(setting==5,"SMA_Type+ Scoliosis+ Age_Enrollment+Y0","SMA_Type+ SMN2_Copy_Number+ Scoliosis+ Age_Enrollment+Y0")
  # subset original data by random indices
  df_b=data[indices,]

  #### Below is the same code as the main function EC_IPW_OPT, except the sandwich SE part
  n=sum(df_b$S)
  m=sum(1-df_b$S)
  pi.S=n/(n+m)

  if(method=="within trial"){

    # estimate ATE
    temp=df_b%>%
      filter(S==1)%>%
      mutate(`piA`=sum(A)/n)%>%
      mutate(w11=`piA`,w10=1-`piA`)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y1","Y2")])

    potential= (temp$A/temp$w11+(1-temp$A)/temp$w10)*Ys
    mu1=colSums(potential[temp$A==1,])/n
    mu0=colSums(potential[temp$A==0,])/n

    tau=mu1-mu0


  }else if(method=="EC-IPW-OPT"){

    piS.model=glm(as.formula(paste("S","~",form_x)) , data =df_b, family = "binomial")

    # estimate ATE
    temp=df_b%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model,newdata =df_b,type="response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`))%>%
      mutate(w11=`piA`,w10=1-`piA`,w00=rx)

    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y1","Y2")])

    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1=colSums(potential[temp$S==1&temp$A==1,])/n
    mu10=colSums(potential[temp$S==1&temp$A==0,])/n
    mu00=colSums(potential[temp$S==0,])/(sum((1-temp$S)*temp$w00))
    mu0=(1-wt)*mu10+wt*mu00

    tau=mu1-mu0


  }else if(method=="EC-AIPW-OPT"){
    piS.model=glm(as.formula(paste("S","~",form_x)) , data =df_b, family = "binomial")
    Y0.model=list(lm(as.formula(paste("Y1","~",form_x)), data = filter(df_b,A==0) ),lm(as.formula(paste("Y2","~",form_x)), data = filter(df_b,A==0) ))


    #for(w in seq(0.05,0.95,by=0.05)){
    # estimate ATE
    temp=df_b%>%
      mutate(`piA`=sum(A[S==1])/n,
             `piS`=sum(S)/(n+m),
             `piSX`=predict(piS.model,newdata =df_b,type="response"),
             rx=(`piSX`/(1-`piSX`))*((1-`piS`)/`piS`),
             `Y1_0`=predict(Y0.model[[1]],newdata =df_b),
             `Y2_0`=predict(Y0.model[[2]],newdata =df_b),
             Y1=Y1-`Y1_0`,
             Y2=Y2-`Y2_0`)%>%
      mutate(w11=`piA`,w10=1-`piA`,w00=rx)

    # check distributional balance:
    #ggplot(dat=data.frame(x=c(rep(1,n),rep(0,m)),y=c(temp$Age_Enrollment[temp$S==1],temp$Age_Enrollment[temp$S==0]),wt=c(rep(1,n),temp$rx[temp$S==0])), aes(y)) + geom_density(aes(weight=wt,group=x))
    #ggplot(dat=data.frame(x=c(rep(1,n),rep(0,m)),y=c(temp$Age_Enrollment[temp$S==1],temp$Age_Enrollment[temp$S==0]),wt=c(rep(1,n),temp$rx[temp$S==0])), aes(y)) + geom_density(aes(group=x))


    ### create outcomes: obs * T
    Ys=as.matrix(temp[c("Y1","Y2")])


    potential= (temp$S*temp$A/temp$w11+temp$S*(1-temp$A)/temp$w10+(1-temp$S)*temp$w00)*Ys
    mu1=colSums(potential[temp$S==1&temp$A==1,])/n
    mu10=colSums(potential[temp$S==1&temp$A==0,])/n
    mu00=colSums(potential[temp$S==0,])/(sum((1-temp$S)*temp$w00))
    mu0=(1-wt)*mu10+wt*mu00

    tau=mu1-mu0



  }


  return(tau)
}


