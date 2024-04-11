#' Simulation toy example with simple covariate strucutrue
#'
#' This function simulate X(baseline covariates), S(trial participation), A(treatment assignemnt), Y(repeated outcomes), in the most basic way, and return a data frame.
#'
#' @param size Sample size 
#' @param tau True treatment effect  
#' @param beta.X Assoiation between X and S
#' @param beta.U Association between U and S, U and Y
#' @return a data frame
#' @examples
#' mysimXY=toy.simXY()
#'  
toy.simXY=function(size=5000, tau=3, beta.X=1, beta.U=-3){
 
  X=rnorm(size, mean=0,sd=2)
  U=rnorm(size, mean=0,sd=2)
  
  piS=exp(log(2)+beta.X*X+beta.U*U)/(1+exp(log(2)+beta.X*X+beta.U*U))
  S=rbinom(size, size=1, prob=piS)
  A=S*rbinom(size, size=1, prob=2/3)
  
  Y=A*tau+X+beta.U*U+rnorm(size,mean = 0,sd=0.5)
  
  simXY=data.frame(X=X,U=U, S=S,A=A,Y=Y)
  return(simXY)
}
