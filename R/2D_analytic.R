#' The 2d analytic density WCP prior for mean and standard deviation of Gaussian distribution
#'
#' @param mean Mean parameter.
#' @param std Standard deviation parameter. 
#' @param eta User specified parameter of the WCP prior.
#' @param base_mean Base model value for the mean parameter.
#'
#' @return A value of density evaluated at c(mean, std).
#' @export
#'
WCP2_2D_Gaussian_analytic = function(mean, 
                                      std,
                                      eta,
                                      base_mean){
  
  if (std < 0 ) {
    stop("std should be a non-negative value!")
  }
  if (mean == Inf | mean == -Inf) {
    stop("mean should be a finite value!")
  }
  density = eta/sqrt((mean - base_mean)^2 + std^2)*exp(-eta * sqrt((mean - base_mean)^2 + std^2))/pi
  return (density)
  
}

#' The 2d analytic density WCP prior for sigma and xi of generalized Pareto distribution
#'
#' @param sigma Sigma parameter.
#' @param xi Xi parameter. 
#' @param eta User specified parameter of the WCP prior.
#'
#' @return A value of density evaluated at c(sigma, xi).
#' @export
#'
WCP1_2D_GP_analytic = function(sigma, 
                                 xi,
                                eta){
  if (sigma < 0 ) {
    stop("sigma should be a non-negative value!")
  }
  if (xi < 0 | xi >= 1) {
    stop("xi should be a value in [0,1)!")
  }
  density = eta/(1 - xi) * exp(-eta * sigma/(1 - xi) )
  return (density)
  
}
