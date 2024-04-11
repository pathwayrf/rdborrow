# Simulate treatment
#' Title
#'
#' @param X 
#' @param S 
#' @param prob 
#'
#' @return Vector A that indicates the treatment status
#' @export
#'
#' @examples
#' A = simulate_trt_assign(X = SyntheticData %>% select(x1, x2), S = SyntheticData %>% select(S), prob = 1/2)
simulate_trt_assign = function(X, S, prob){
  # sanity check: 
  # TODO: X, S has same number of rows
  n = nrow(X)
  A = S$S * rbinom(n, size=1, prob = prob)
  A = data.frame(A)
}