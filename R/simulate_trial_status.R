# trial participation status
## TODO: generalize to incorporate other models
#' Simulate trial status indicator
#'
#' @param X 
#' @param model_specs 
#'
#' @return a data frame containing the trial status vector
#' @export
#'
#' @examples
#' S = simulate_trial_status(X, model_specs = list(
#'   family = "binomial",
#'   coef = c(1,2,3)
#'   ))
simulate_trial_status = function(X, model_specs){
  n = nrow(X)
  X_intercept = cbind(intercept = 1, X)
  if (model_specs$family == "binomial"){
    beta = model_specs$coef
    piS = inv.logit(as.matrix(X_intercept) %*% beta)
    #    piS = exp(log(2)+beta.X*X+beta.U*U)/(1+exp(log(2)+beta.X*X+beta.U*U))
  }
  S = rbinom(n, size=1, prob=piS)
  data.frame(S = S)
}