#' method class
#'
#' @slot model_form_mu0_ext 
#' @slot model_form_mu0_rct 
#' @slot model_form_mu1_rct 
#' @slot model_form_piS character. 
#' @slot bootstrap_flag 
#' @slot model_form_piA 
#' @slot bootstrap_obj 
#'
#' @include method_class.R
#' @include bootstrap_class.R
#' @export setup_method_DID
#'
#' @examples
#' method_DID_obj = setup_method_DID(
#'    method_name = "IPW",
#'    bootstrap_flag = T,
#'    bootstrap_obj = bootstrap_obj,
#'    model_form_piS = "S ~ x1 + x2 + x3 + x4 + x5",
#'    model_form_piA = "A ~ x1 + x2 + x3 + x4 + x5")
#' 
.method_DID_obj = setClass(
  "method_DID_obj",
  contains = "method_OLE_obj",
  slots = c(
    bootstrap_flag = "logical",
    bootstrap_obj = "bootstrap_obj",
    model_form_piS = "character",
    model_form_piA = "character",
    model_form_mu0_ext = "character",
    model_form_mu0_rct = "character",
    model_form_mu1_rct = "character"
  ),
  prototype = list(
    bootstrap_flag = F,
    bootstrap_obj = .bootstrap_obj(),
    model_form_piA = "",
    model_form_piS = "",
    model_form_mu0_ext = "",
    model_form_mu0_rct = "",
    model_form_mu1_rct = ""
  )
)


setup_method_DID = function(method_name = "IPW", 
                            bootstrap_flag = F,
                            bootstrap_obj = .bootstrap_obj(),
                            model_form_piS = "",
                            model_form_mu0_ext = "",
                            model_form_piA = "",
                            model_form_mu0_rct = "",
                            model_form_mu1_rct = ""){
  # TODO: sanity check
  # correct initialization of objects
  # correct dimension compatible
  # validity
  # model_form_mu0 and the dimension of the random vector
  
  method_DID_obj = .method_DID_obj(
    method_name = method_name,
    bootstrap_flag = bootstrap_flag,
    bootstrap_obj = bootstrap_obj,
    model_form_piS = model_form_piS,
    model_form_mu0_ext = model_form_mu0_ext,
    model_form_piA = model_form_piA,
    model_form_mu0_rct = model_form_mu0_rct,
    model_form_mu1_rct = model_form_mu1_rct
  )
  
  
}