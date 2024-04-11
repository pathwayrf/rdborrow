#' Simulate X from a mixture model
#' @param n total number of units simulated
#' @param p_cat dimension of categorical covariates
#' @param p_cont dimension of continuous covariates
#' @param cat_level_list a list describing the levels of categorical variables
#' @param cat_comb_prob probability of each combination of categorical variables
#' @return a list contains simulated covariates
#' @export
#' 
#' @examples
#' X = simulate_X_mixture(n = 100, p_cat = 0, p_cont = 2, 
#' cat_level_list = list(),
#' cat_comb_prob = c(),
#' cont_para_list = list(list(mean = c(0, 0), sigma = diag(2)))
#' )



simulate_X_mixture = function(n, p_cat, p_cont, cat_level_list, cat_comb_prob, cont_para_list){
  # TODO: sanity check
  # p_cat and cat_level_list
  # cat_level_list and cat_comb_prob
  # cat_comb_prob is a prob vector
  # p_cont and cont_para_list
  
  p = p_cat + p_cont
  
  if (p_cat == 0){
    covariate = data.frame(rmvnorm(n, mean = cont_para_list[[1]]$mean, sigma = cont_para_list[[1]]$sigma))
  }
  else{
    # TODO: extend to general distributions
    print(cbind(expand.grid(cat_level_list), p = cat_comb_prob))
    
    # TODO: generate description of the full distribution
    N = c(rmultinom(1, n, cat_comb_prob))
    num_comb = length(cat_comb_prob)
    
    covariate_cat = expand.grid(cat_level_list)
    covariate_cat = covariate_cat %>% mutate(count = N) %>% uncount(count)
    # flog.debug(paste("User chose to simulate continuous covariate. \n"))
    
    if (p_cont == 0){
      covariate = covariate_cat
    }
    else{
      covariate_cont =  bind_rows(lapply(1:num_comb, function(k){data.frame(rmvnorm(N[k], mean =    cont_para_list[[k]]$mean, sigma = cont_para_list[[k]]$sigma))}))
      
      covariate = cbind(covariate_cat, covariate_cont)
    }
    
  }
  
  
  colnames(covariate) = paste0("x", 1:p)
  covariate
}
