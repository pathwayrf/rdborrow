#' simulate X by coupling several marginal distributions using copula
#' @param n total number of units simulated
#' @param p dimension of the parameters
#' @param cp copula
#' @param margins marginal distributions
#' @param paramMargins parameters for the marginal distributions
#' @return a list contains simulated data and true ATE
#' @export
#' 
#' @examples
#' normal <- normalCopula(param = c(0.8), dim = 4, dispstr = "ar1")
#' X <- simulate_X_copula(1000, 4, normal, 
#'                       margins = c("norm", "t", "norm", "binom"), 
#'                       paramMargins = list(list(mean = 2, sd=3),
#'                                  list(df = 2),
#'                                  list(mean = 0, sd = 1),
#'                                  list(size = 10, prob = 0.5))
#'                       )
#' cor(X, method = "spearman")


simulate_X_copula = function(n, p, cp, margins, paramMargins){
  # TODO: sanity check
  # p and cp
  
  multivariate_dist <- mvdc(copula = cp,
                            margins = margins,
                            paramMargins = paramMargins)
  
  
  # flog.debug(paste("trial has sample size =", n, "\n"))
  # print(multivariate_dist)
  
  covariate = rMvdc(n, multivariate_dist)
  # flog.debug(paste("User chose to simulate continuous covariate. \n"))
  
  colnames(covariate) = paste0("x", 1:p)
  data.frame(covariate)
}