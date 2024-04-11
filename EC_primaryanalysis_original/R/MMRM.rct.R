#' Run MMRM on RCT data
#'
#' This function runs the MMRM model on the RCT data and returns estimate, SE and confidence intervals for ATE at first two time points
#'
#' @param setting A numeric value.
#' Setting 1: all confounder correctly specified
#' Setting 2: outcome model misspecifed
#' Setting 3: propensity model misspecifed
#' Setting 4: both propensity model and outcome models misspecified
#' Setting 5: Exists unmeasured confounder: SMM2_Copy number as pseudo omitted confounder
#' @param df A data frame contains simulated data, returned by cal.simXY()
#' @return A list contains estimate of ATE (point estimate, SE, lower, upper)
#' @examples
#' result=cal.simXY()
#' list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
#' MMRM.rct(df=mysimXY)



MMRM.rct<-function(setting=1,
                      df){
  df_l=df%>%
    tidyr::gather(key="var",value="y",Y1:Y2 )%>%
    mutate(time=substr(var,2,2))%>%
    arrange(id,time)%>%
    dplyr::filter(S==1)%>%
    mutate(time=factor(time),
           id=as.factor(id))



  library(mmrm)
  if(setting!=5){
    m1 <- mmrm(
      formula = y ~ time + A+ time*A +SMA_Type+ SMN2_Copy_Number+Scoliosis+Age_Enrollment +Y0+ us(time | id),
      data = df_l
    )
    summary(m1)

    contrast <- numeric(length(component(m1, "beta_est")))
    contrast[3] <- 1
    tau1=df_1d(m1, contrast)

    contrast <- numeric(length(component(m1, "beta_est")))
    contrast[3]=contrast[9] <- 1
    tau2=df_1d(m1, contrast)

  }else{
    m1 <- mmrm(
      formula = y ~ time + A+ time*A +SMA_Type+Scoliosis+Age_Enrollment +Y0+ us(time | id),
      data = df_l
    )
    summary(m1)

    contrast <- numeric(length(component(m1, "beta_est")))
    contrast[3] <- 1
    tau1=df_1d(m1, contrast)

    contrast <- numeric(length(component(m1, "beta_est")))
    contrast[3]=contrast[8] <- 1
    tau2=df_1d(m1, contrast)

  }


  return(list(c(tau1$est,tau1$se,tau1$est-1.96*tau1$se,tau1$est+1.96*tau1$se),
              c(tau2$est,tau2$se,tau2$est-1.96*tau2$se,tau2$est+1.96*tau2$se)))

}


# result=cal.simXY()
# list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
# MMRM.rct(df=mysimXY)


