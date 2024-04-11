#' Simulate outcome from given effect additive models
#'
#' @param n 
#' @param T_follow 
#' @param X 
#' @param A 
#' @param outcome_model_specs 
#' 
#' @return a data frame containing simulated outcome 
#' @export
#'
#' @examples 
#' 
#' model_form1 = "y1 = A*3 + x1*1 + x2*1 + rnorm(n, mean = 0, sd=0.5)"
#' model_form2 = "y2 = A*0 + x1*1 + x2*(-1) + rnorm(n, mean = 0, sd=0.5)"
#' Y = simulate_outcome_from_model(T_follow = 2, X = X, A = A, 
#'                                 outcome_model_specs = list(
#'                                   list(
#'                                     model_form = model_form1
#'                                   ),
#'                                   list(
#'                                     model_form = model_form2
#'                                   )
#'                                 ))
#' Y

simulate_outcome_from_model = function(X, A, outcome_model_specs, OLE_flag, T_cross){
  # TODO: sanity check 
  # T_cross > 0
  # OLE_flag and length of T_cross
  n = nrow(X)
  T_follow = length(outcome_model_specs)
  
  if (!OLE_flag){
    full_data = data.frame(X, A)
    num_col = ncol(full_data)
    
    for (t in 1:T_follow){
      # print("hi")
      Y = with(full_data, {
        with(outcome_model_specs[[t]], {
          # print("hihi")
          model_form_full = paste0("y = A*", effect, " + " , model_form_x)
          # print(model_form_full)
          Y = eval(parse(text = model_form_full))
          Y + rnorm(n, mean = noise_mean, sd = noise_sd)
        })
      })
      full_data = cbind(full_data, Y)
      colnames(full_data)[ncol(full_data)] = paste0("y", t)
    }
  }else{
    full_data = data.frame(X, A)
    num_col = ncol(full_data)
    for (t in 1:T_cross){
      Y = with(full_data, {
        with(outcome_model_specs[[t]], {
          model_form_full = paste0("y = A*", effect, " + " , model_form_x)
          Y = eval(parse(text = model_form_full))
          Y + rnorm(n, mean = noise_mean, sd = noise_sd)
        })
      })
      full_data = cbind(full_data, Y)
      colnames(full_data)[ncol(full_data)] = paste0("y", t)
      attr(full_data[, ncol(full_data)], "label") = paste0("Period ", t)
    }
    
    for (t in (T_cross+1):T_follow){
      Y = with(full_data, {
        with(outcome_model_specs[[t]], {
          model_form_full = paste0("y = 1*", effect, " + " , model_form_x)
          Y = eval(parse(text = model_form_full))
          Y + rnorm(n, mean = noise_mean, sd = noise_sd)
        })
      })
      full_data = cbind(full_data, Y)
      colnames(full_data)[ncol(full_data)] = paste0("y", t)
      attr(full_data[, ncol(full_data)], "label") = paste0("Period ", t)
    }
    
  }
  
  
  
  # return outcome data from simulation
  full_data[(num_col + 1):(num_col + T_follow)]
  
}