#' Implement the (naive) Reference-based multiple imputation 
#'
#' 
#' @param df A data frame contains simulated data, returned by cal.simXY()
#' @param setting The setting number
#' Setting 1: No unmeasured confounder 
#' Setting 2: unmeasured confounding with time-invariant effect, all models with observed covariates are correctly specified
#' Setting 3: unmeasured confounding with time-invariant effect, outcome models with observed covariates are correctly specified
#' Setting 4: unmeasured confounding with time-invariant effect, PS models with observed covariates are correctly specified
#' Setting 5: unmeasured confounding with time-invariant effect, both models with observed covariates are slightly misspecified
#' Setting 6: unmeasured confounding with time-variant effect, both models with observed covariates are slightly misspecified
#' @return A list contains: estimated ATE, SE, weihgt used, SE by Bootstrap and a 95% confidence interval for primary endpoint (only when Bootstrap=TRUE)
#' @examples
#' RBMI(df=mysimXY,setting=1)
#' 
 

RBMI = function(df,
             setting){
  
  # covfill: all variables available about each subject: X and S, A and Study, 
  # covim: covariate for imputation model
  # covmod: covariate included in the final ancova analysis
  
  if(setting==1){
    covfill=c("SMA_Type", "Scoliosis", "SMN2_Copy_Number", "Age_Enrollment", "Y0","Y1","Y2","S","A","study")
    covim=c("SMA_Type*visit", "Scoliosis*visit", "SMN2_Copy_Number*visit", "Age_Enrollment*visit", "Y0*visit")
    
    covmod=c("SMA_Type", "Scoliosis", "SMN2_Copy_Number", "Age_Enrollment", "Y0")
  }else{
    covfill=c("SMA_Type", "Scoliosis","SMN2_Copy_Number", "Age_Enrollment", "Y0","Y1","Y2","S","A","study")
    covim=c("SMN2_Copy_Number*visit", "Scoliosis*visit", "Age_Enrollment*visit", "Y0*visit")
    covmod=c("SMN2_Copy_Number", "Scoliosis", "Age_Enrollment", "Y0")
  }
  
  dat <- df%>% 
    mutate(S=factor(S))%>%
    mutate(A=factor(A))%>%
    mutate(study=factor(study))%>%
    pivot_longer(c(Y1,Y2,Y3,Y4), names_to = "key", values_to = "Y")%>%
    mutate(visit=factor(substr(key,2,2)))%>%
    filter(!(study=="RCT Control Arm" & visit %in%c(3,4)))%>%
    dplyr::select(-"key")
  
  dat<-dat%>%left_join(dplyr::select(df,id,Y1,Y2),by=c("id"))%>%
    mutate(id=factor(id))
  
  
  # Use expand_locf to add rows corresponding to visits with missing outcomes to the dataset
  dat <- expand_locf(
    dat,
    id = levels(dat$id), # expand by PATIENT and VISIT 
    visit= levels(dat$visit),
    vars = c(covfill), # fill with LOCF BASVAL and THERAPY
    group = c("id"),
    order = c("id", "visit")
  )
  #View(dat)
  
  
  
  ################ Draws
  
  # create data_ice and set the imputation strategy to JR for
  # each patient with at least one missing observation
  dat_ice <- dat%>% 
    arrange(id,visit) %>% 
    filter(is.na(Y)) %>% 
    group_by(id) %>% 
    dplyr::slice(1)%>%
    ungroup() %>% 
    dplyr::select(id,visit) %>% 
    mutate(strategy = "JR")
  
  
  #%>%mutate(visit=factor(visit))
  dat <- expand_locf(
    dat,
    id = levels(dat$id), # expand by PATIENT and VISIT 
    visit= levels(dat$visit),
    vars = c(covfill), # fill with LOCF BASVAL and THERAPY
    group = c("id"),
    order = c("id", "visit")
  )
  
  
  # Define the names of key variables in our dataset and
  # the covariates included in the imputation model using `set_vars()`
  # Note that the covariates argument can also include interaction terms
  vars <- set_vars(
    outcome = "Y",
    visit = "visit",
    subjid = "id",
    group = "study",
    covariates = covim
  )
  
  # Define which imputation method to use (here: Bayesian multiple imputation with 150 imputed datsets) 
  method <- method_bayes(
    burn_in = 200,
    burn_between = 5,
    n_samples = 500,
    seed = 675442751
  )
  
  
  # Create samples for the imputation parameters by running the draws() function
  #set.seed(987)
  drawObj <- draws(
    data = dat,
    data_ice = dat_ice,
    vars = vars,
    method = method,
    quiet = TRUE
  )
  drawObj
  
  ################ Impute
  
  imputeObj <- impute(
    drawObj,
    references = c("RCT Control Arm"= "External Control")
  )
  imputeObj 
  
  imputed_dfs <- extract_imputed_dfs(imputeObj)
  
  #head(imputed_dfs[[10]], 12) # first 12 rows of 10th imputed dataset
  # remove external control
  # rct.trt <- df%>% 
  #   filter(A==1)%>%
  #   mutate(S=factor(S))%>%
  #   mutate(A=factor(A))%>%
  #   pivot_longer(c(Y1,Y2,Y3,Y4), names_to = "key", values_to = "Y")%>%
  #   mutate(visit=factor(substr(key,2,2)))%>%
  #   filter(!(study=="RCT Control Arm" & visit %in%c(3,4)))%>%
  #   dplyr::select(-"key")%>%left_join(dplyr::select(df,id,Y1,Y2),by=c("id"))%>%
  #   mutate(id=factor(id))
  
  for(i in 1:length(imputed_dfs)){imputed_dfs[[i]]=filter(imputed_dfs[[i]],S==1)}
  
  ############# Analyse
  
  anaObj <- analyse(
    imputeObj,
    ancova,
    vars = set_vars(
      subjid = "id",
      outcome = "Y",
      visit = "visit",
      group = "A",
      covariates = covmod
    )
  )
  anaObj
  
  ################ Pool
  poolObj <- pool(
    anaObj, 
    conf.level = 0.95, 
    alternative = "two.sided"
  )
  return(list(c(poolObj$pars$trt_3$est, poolObj$pars$trt_3$se, poolObj$pars$trt_3$ci),
              c(poolObj$pars$trt_4$est, poolObj$pars$trt_4$se, poolObj$pars$trt_4$ci)))
  
}
 



