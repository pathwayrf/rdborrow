#' method class
#'
#' @slot method_name character. 
#' @slot bootstrap_flag 
#' @slot bootstrap_obj 
#'
#' @include bootstrap_class.R
#' @export setup_method
#'
#' @examples
#' method_weighting_obj = setup_method_weighting(
#'   method_name = "AIPW",
#'   optimal_weight_flag = T, 
#'   wt = 0,
#'   bootstrap_flag = T,
#'   bootstrap_obj = bootstrap_obj,
#'   model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5",
#'   model_form_mu0_ext = c("y1 ~ x1 + x2 + x3 + x4 + x5",
#'                          "y2 ~ x1 + x2 + x3 + x4 + x5"))
#' 

.method_obj = setClass(
  "method_obj",
  slots = c(
    method_name = "character", # method_name
    bootstrap_flag = "logical",
    bootstrap_obj = "bootstrap_obj"
  )
)

.method_primary_obj = setClass(
  "method_primary_obj",
  contains = "method_obj"
)

.method_OLE_obj = setClass(
  "method_OLE_obj",
  contains = "method_obj"
)


setup_method = function(method_name = "", 
                        bootstrap_flag = F,
                        bootstrap_obj = .bootstrap_obj()){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  # model_form_mu0 and the dimension of the random vector
  
  method_obj = .method_obj(
    method_name = method_name,
    bootstrap_flag = bootstrap_flag,
    bootstrap_obj = bootstrap_obj
  )
  
  
}