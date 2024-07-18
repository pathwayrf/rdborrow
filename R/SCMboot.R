SCMboot = function(data, 
                   indices,
                   outcome_col_name, 
                   trial_status_col_name, 
                   treatment_col_name, 
                   covariates_col_name,
                   T_cross,
                   lambda,
                   pb){
  # browser()
  pb$tick()
  # print(Sys.getpid())
  # plan(multicore)
  df_b = data[indices, ]
  
  n10 = sum(df_b$S==1 & df_b$A==0) # number of rct control subjects
  n00 = sum(df_b$S==0) # number of external control subjects
  
  Y = subset(df_b, select = outcome_col_name)
  S = subset(df_b, select = trial_status_col_name)
  A = subset(df_b, select = treatment_col_name)
  X = subset(df_b, select = covariates_col_name)
  names(S) = "S"
  names(A) = "A"
  
  # standardize numeric vars
  # columns_to_standardize = c("SMN2_Copy_Number", "Age_Enrollment", "Y0")
  # data[columns_to_standardize] = lapply(data[columns_to_standardize], scale)
  
  # extract long term column names
  short_term_col_name = outcome_col_name[1:T_cross]
  long_term_col_name = outcome_col_name[(T_cross + 1):length(outcome_col_name)]
  
  # create data matrices: attributes by row and subject by column
  X10 = t(as.matrix(df_b %>% filter(S==1 & A==0) %>% dplyr::select(all_of(c(covariates_col_name, outcome_col_name)))))
  X00 = t(as.matrix(df_b %>% filter(S==0) %>% dplyr::select(all_of(c(covariates_col_name, outcome_col_name)))))
  
  # remove colnames of X10 and X00
  colnames(X00) = NULL
  colnames(X10) = NULL
  
  # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
  res = lapply(1:n10, subject_SC, X10 = X10, X00 = X00, long_term_col_name = long_term_col_name, lambda = lambda)
  
  # weight matrix: n0*m (combinging by rows)
  # wt.mat = do.call(rbind, lapply(res, function(x) as.vector(x[[1]])))
  
  # estimated outcome matrix: n0*T (combining by rows)
  y.est.mat = do.call(rbind, lapply(res, function(x) as.vector(x[[2]])))
  
  # Aggregate group level synthetic control estiamte
  Y.trt = df_b %>% filter(S == 1 & A == 1) %>% 
    dplyr::select(all_of(long_term_col_name)) %>% colMeans()
  
  tau = Y.trt - colMeans(y.est.mat)
  
  return(tau)
}


