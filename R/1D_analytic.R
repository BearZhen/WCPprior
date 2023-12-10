#' The 1d analytic density WCP prior for phi of stationary AR1 process
#'
#' @param phi phi parameter.
#' @param eta User specified parameter of the WCP prior.
#' @param n Length of AR1 process. 
#' @param sigma Standard deviation of the process.
#'
#' @return A value of density evaluated at c(mean, std).
#' @export
#'
WCP_1D_AR1_analytic = function(phi, 
                               eta,
                               n,
                               sigma){
  if (n < 1 ) {
    stop("n should be no less than 1!")
  }
  if (phi < -1 | phi >=1 ) {
    stop("phi should be in [-1, 1)!")
  }
  if (sigma <= 0 ) {
    stop("sigma should be positive!")
  }
  
  f = sqrt(n*(1 - phi^2) - 2*phi*(1 - phi^n))
  c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) )
  J = sigma/sqrt(2) * ( (n*phi^n - 1 + phi^n - n*phi)*(1 - phi) + f^2 )/( f*sqrt(n - f/(1-phi))*(1 - phi)^2 )
  density = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - phi) )  )/(1 - exp(-eta*c))
  return (density)
  
}


#' The 1d analytic density WCP prior for xi of generalized Pareto distribution
#'
#' @param xi Xi parameter. 
#' @param eta User specified parameter of the WCP prior.
#'
#' @return A value of density evaluated at c(sigma, xi).
#' @export
#'
WCP_2D_GP_analytic = function(xi,
                              eta){
  if (xi < 0 | xi >= 1) {
    stop("xi should be a value in [0,1)!")
  }
  density = eta/(1 - xi)^2 * exp( - eta * xi/(1 - xi) )
  return (density)
}

  