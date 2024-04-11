#' simulate X by discretizing a multivariate normal distribution
#' pkg: futile.logger, mvtnorm
#' @param n total number of units simulated
#' @param p dimension of the parameters
#' @param mu mean for mvnorm
#' @param sig cov of mvnorm
#' @param cat_cols columns with categorical variables
#' @param cat_prob list
#' 
#' @return a list contains simulated data and true ATE
#' @export
#' 
#' @examples
#' simulate_X_dct_mvnorm(20, 3, mu = rep(0, 3), sig = diag(3), cat_cols = c(1), cat_prob = list(c(0.3, 0.7)))
#' simulate_X_dct_mvnorm(20, 3, mu = rep(0, 3), sig = diag(3), cat_cols = c(1, 3), cat_prob = list(c(0.2, 0.6, 0.2), c(0.3, 0.7)))


simulate_X_dct_mvnorm = function(n, p, mu = rep(0, p), sig = diag(p), cat_cols = c(), cat_prob = list()){
  # TODO: sanity check
  # cat_col: numeric vector
  # cat_prob: a list, same length as cat_col, each element has sum 1 and > 0
  p_cat = length(cat_cols)
  flog.debug(paste("trial has sample size =", n, "\n"))
  # cov_st = sum(grepl("cov", colnames(dt))) + 1
  
  covariate = NULL
  # flog.debug(paste("[s_cov] Total number of covariates to add this time =", n_cov, "\n"))
  
  covariate = rmvnorm(n, mean = mu, sigma = sig)
  # flog.debug(paste("User chose to simulate continuous covariate. \n"))
  if (p_cat > 0) {
    covariate[, cat_cols] = sapply(1:p_cat, function(k){
      cut_val = qnorm(cumsum(cat_prob[[k]]), mean = mu[cat_cols[k]], sd = sqrt(sig[cat_cols[k], cat_cols[k]]))
      # flog.debug(cat("User chose to simulate", n_cov, "correlated binary covariate(s) with cut_val =", cut_val, "\n"))
      
      sapply(1:n, function(i){ which(covariate[i, cat_cols[k]] < cut_val)[1] - 1 })
    })
  }
  
  
  colnames(covariate) = paste0("x", 1:p)
  data.frame(covariate)
}