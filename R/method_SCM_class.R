#' method class
#'
#' @slot parallel 
#' @slot ncpus 
#' @slot lambda.min 
#' @slot lambda.max 
#' @slot nlambda 
#'
#' @include method_class.R
#' @include bootstrap_class.R
#' @export setup_method_SCM
#'
#' @examples
#' method_SCM_obj = setup_method_SCM(
#'   method_name = "SCM", 
#'   bootstrap_flag = T, 
#'   bootstrap_obj = bootstrap_obj, 
#'   lambda.min = 0, 
#'   lambda.max = 1e-3, 
#'   nlambda = 10, 
#'   parallel = "multicore", 
#'   ncpus = 4) 
#' 

.method_SCM_obj = setClass(
  "method_SCM_obj",
  contains = "method_OLE_obj",
  slots = c(
    lambda.min = "numeric", # minimum value 
    lambda.max = "numeric", 
    nlambda = "numeric", # nfolds
    parallel = "character", 
    ncpus = "numeric"
  ),
  prototype = list(
    nlambda = 10, 
    parallel = "no", 
    ncpus = 1
  )
)


setup_method_SCM = function(method_name = "SCM", 
                            bootstrap_flag = F,
                            bootstrap_obj = .bootstrap_obj(),
                            lambda.min,
                            lambda.max,
                            nlambda = 10, 
                            parallel = "no", 
                            ncpus = 1){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  # model_form_mu0 and the dimension of the random vector
  
  method_SCM_obj = .method_SCM_obj(
    method_name = method_name, 
    bootstrap_flag = bootstrap_flag,
    bootstrap_obj = bootstrap_obj,
    lambda.min = lambda.min,
    lambda.max = lambda.max,
    nlambda = nlambda, 
    parallel = parallel, 
    ncpus = ncpus
  )
  
  
}