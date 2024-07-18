#' Implement the Synthetic Control Method  
#'
#' SCM() is the main function calculates the estimated ATE by SC method and Bootstrap CI, it calls subject_SC() and lambdacv().
#'
#' @param outcome_col_name 
#' @param trial_status_col_name 
#' @param treatment_col_name 
#' @param covariates_col_name 
#' @param T_cross 
#' @param Bootstrap 
#' @param R 
#' @param bootstrap_CI_type 
#' @param alpha 
#' @param data A data frame
#' @param lambda.min 
#' @param lambda.max 
#' @param nlambda 
#' @param parallel 
#' @param ncpus 
#'
#' @include SCMboot.R
#' @return A list contains: estimated ATE, SE, weight used, SE by Bootstrap and a 95% confidence interval for primary endpoint (only when Bootstrap=TRUE)
#' @examples This is part of the run_analysis() function
#' 
#'  
SCM = function(data,
               outcome_col_name, 
               trial_status_col_name, 
               treatment_col_name, 
               covariates_col_name,
               T_cross,
               Bootstrap = T,
               R = 1e2,
               bootstrap_CI_type = "bca", 
               alpha = 0.05,
               lambda.min = 0,
               lambda.max = 0.1,
               nlambda = 2,
               parallel = "no",
               ncpus = 1,
               quiet = TRUE){
  
  cat("Running the synthetic control method...\n")
  
  Y = subset(data, select = outcome_col_name)
  S = subset(data, select = trial_status_col_name)
  A = subset(data, select = treatment_col_name)
  X = subset(data, select = covariates_col_name)
  names(S) = "S"
  names(A) = "A"
  
  df = data.frame(Y, S = S, A = A, X)
  
  n10 = sum(df$S==1 & df$A==0) # number of rct control subjects
  n00 = sum(df$S==0) # number of external control subjects
  T_follow = length(outcome_col_name)
  
  # standardize numeric vars
  # columns_to_standardize = c("SMN2_Copy_Number", "Age_Enrollment", "Y0")
  # df[columns_to_standardize] = lapply(df[columns_to_standardize], scale)
  
  # extract long term column names
  short_term_col_name = outcome_col_name[1:T_cross]
  long_term_col_name = outcome_col_name[(T_cross + 1):length(outcome_col_name)]
  
  # create df matrices: attributes by row and subject by column
  X10 = t(as.matrix(df %>% filter(S==1 & A==0) %>% dplyr::select(all_of(c(covariates_col_name, outcome_col_name)))))
  X00 = t(as.matrix(df %>% filter(S==0) %>% dplyr::select(all_of(c(covariates_col_name, outcome_col_name)))))
  
  # remove colnames of X10 and X00
  colnames(X00) = NULL
  colnames(X10) = NULL
  
  
  # Use LOOCV to find the optimal penalty parameter lambda
  cat("Performing cross validation for tuning parameter selection...\n")
  
  pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                        total = nlambda, 
                        clear = FALSE,
                        show_after = 0)
  
  pb$tick(len = 0)   # now it displays the 0% bar
  
  lambda = lambdacv(ec = X00,
                    long_term_col_name = long_term_col_name,
                    lambda.min = lambda.min,
                    lambda.max = lambda.max,
                    nlambda = nlambda,
                    pb = pb)
  
  cat("Constructing pseudo controls for internal data...\n")
  # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
  res = lapply(1:n10, 
               subject_SC, 
               X10 = X10, 
               X00 = X00, 
               long_term_col_name = long_term_col_name, 
               lambda = lambda)
  
  # weight matrix: n0*m (combinging by rows)
  wt.mat = do.call(rbind, lapply(res, function(x) as.vector(x[[1]])))
  
  # estimated outcome matrix: n0*T (combining by rows)
  y.est.mat = do.call(rbind, lapply(res, function(x) as.vector(x[[2]])))
  
  # Aggregate group level synthetic control estiamte
  Y.trt = df %>% filter(S == 1 & A == 1) %>% 
    dplyr::select(all_of(long_term_col_name)) %>% colMeans()

  tau = Y.trt - colMeans(y.est.mat)
  names(tau) = paste0("tau", (T_cross+1):T_follow)
  
  ####### Use Bootstrap for standard error and confidence intervals
  cat("Performing bootstrap inference with SCM estimates...\n")
  
  Group_ID = df %>% group_by(S, A) %>% mutate(group_id = cur_group_id())
  Group_ID = Group_ID$group_id
  
  boot.ci.type = switch(bootstrap_CI_type,
                         norm = "normal",
                         bca = "bca",
                         stud = "student",
                         perc = "percent",
                         basic = "basic"
  )
  
  pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                        total = R + 1,
                        clear = FALSE,
                        show_after = 0)
  pb$tick(len = 0)   # now it displays the 0% bar
  
  start_time = Sys.time()
  boot.out = boot(data = df,
                   statistic = SCMboot,
                   outcome_col_name = outcome_col_name, 
                   trial_status_col_name = trial_status_col_name, 
                   treatment_col_name = treatment_col_name, 
                   covariates_col_name = covariates_col_name,
                   T_cross = T_cross,
                   lambda = lambda,
                   pb = pb,
                   parallel = parallel,
                   ncpus = ncpus,
                   R = R,
                   strata = Group_ID)
  end_time = Sys.time()
  cat("time elapsed for bootstrap: ", round(end_time - start_time, 2))
  
  lower_CI_boot = sapply(1:(T_follow - T_cross), 
                         function(x) {
                           ci = boot.ci(boot.out, conf = 1 - alpha, type = bootstrap_CI_type, index = x)
                           ci[[boot.ci.type]][4]
                         })
  
  upper_CI_boot = sapply(1:(T_follow - T_cross), 
                         function(x) {
                           ci = boot.ci(boot.out, conf = 1 - alpha, type = bootstrap_CI_type, index = x)
                           ci[[boot.ci.type]][5]
                         })
  
  sd.boot = sqrt(diag(var(boot.out$t)))
  
  results = data.frame(
    point_estimates = tau,
    lower_CI_boot = lower_CI_boot,
    upper_CI_boot = upper_CI_boot
  )
  
  return(results)
  
}


#' find synthetic control for one specific subject
#'
#' @param subject 
#' @param X10 
#' @param X00 
#' @param long_term_col_name 
#' @param lambda 

subject_SC = function(subject, X10, X00, long_term_col_name, lambda){
  # subset to only the intersted subject and matching vars
  x1 = X10[-which(row.names(X10) %in% long_term_col_name), subject]
  X0 = X00[-which(row.names(X10) %in% long_term_col_name), ]
  
  # optimization
  w = Variable(dim(X00)[2])
  loss = sum(((x1 - X0%*%w))^2)
  penal = function(w, x1, X0, lambda = lambda){
    lambda * (sum(colSums((x1-X0)^2)*w))
  }
  obj = loss + penal(w = w, x1 = x1, X0 = X0, lambda = lambda)
  constr = list(sum(w) == 1, w >= 0)
  prob = Problem(Minimize(obj), constr)
  result = solve(prob, solver = 'ECOS') # OSQP, SCS, ECOS
  # estimated weight for this subject
  wt.est = result$getValue(w)
  # sc estimate for all time points for this subject
  y.est = X00[long_term_col_name, ] %*% wt.est
  
  return(list(wt.est, y.est))
}

#' find the optima lambda via LOOCV 
#'
#' @param ec 
#' @param lambda.min 
#' @param lambda.max 
#' @param nlambda 
#' @param long_term_col_name 

lambdacv = function(ec,
                    long_term_col_name,
                    lambda.min = 0,
                    lambda.max = 0.1,
                    nlambda = 10,
                    pb = NULL){
  # print(pb)
  pb$tick(0)
  lambda_vals = seq(lambda.min, lambda.max, length.out = nlambda)
  # leave-one-out cv error for individul external controls
  mse_vals = sapply(lambda_vals,
                    function (lambda) {
                      pb$tick()
                      res = lapply(1:dim(ec)[2], function(loocv){
                        # subset to only the intersted subject and matching vars
                        x1 = ec[-which(row.names(ec) %in% long_term_col_name), loocv]
                        X0 = ec[-which(row.names(ec) %in% long_term_col_name), -loocv]
                        
                        # optimization
                        w = Variable(dim(ec)[2]-1)
                        loss = sum(((x1 - X0%*%w))^2)
                        penal = function(w, x1, X0, lambda = lambda){lambda*(sum(colSums((x1-X0)^2)*w))}
                        obj = loss + penal(w=w, x1=x1, X0=X0, lambda = lambda)
                        constr = list(sum(w) == 1, w>= 0)
                        prob = Problem(Minimize(obj), constr)
                        result = solve(prob, solver = 'ECOS') # OSQP, SCS, ECOS
                        # estimated weight for this subject
                        wt.est = result$getValue(w)
                        # sc estimate for all time points for this subject
                        y.est = ec[long_term_col_name, -loocv] %*% wt.est
                        
                        return(list(wt.est, y.est))
                      })
                      # estimated outcome matrix: n0*T (combining by rows)
                      y.est.mat = do.call(rbind, lapply(res, function(x) as.vector(x[[2]])))
                      
                      # calcualte MSE
                      mean((ec[long_term_col_name, ] - t(y.est.mat))^2)
                      
                    })
  
  return(lambda_vals[which.min(mse_vals)])
  
}

