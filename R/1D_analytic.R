#' The 1d analytic density WCP prior for phi of stationary AR1 process
#'
#' @param seq_phi A vector of values of phi parameter.
#' @param eta User specified parameter of the WCP prior.
#' @param n Length of AR1 process. 
#' @param sigma Standard deviation of the process.
#' @param inla_table Should the results be returned as a table to be used in INLA? Default is FALSE.
#'
#' @return A list including densities on seq_phi and a inla table
#' @export
#'
WCP2_1D_AR1_analytic = function(seq_phi, 
                               eta,
                               n,
                               sigma,
                               inla_table = FALSE){
  if (n < 1 ) {
    stop("n should be no less than 1!")
  }
  if (sigma <= 0 ) {
    stop("sigma should be positive!")
  }
  density = numeric()
  for (i in 1:length(seq_phi)){
    phi = seq_phi[i]
    if (phi < -1 | phi > 1 ) {
      stop("phi should be in [-1, 1]!")
    }
    f = sqrt(n*(1 - phi^2) - 2*phi*(1 - phi^n))
    c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) )
    J = sigma/sqrt(2) * ( (n*phi^n - 1 + phi^n - n*phi)*(1 - phi) + f^2 )/( f*sqrt(n - f/(1-phi))*(1 - phi)^2 )
    density[i] = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - phi) )  )/(1 - exp(-eta*c))
  }
  if (inla_table){
  # inla table
  inla.prior.table = paste0("table: ",
                        paste(c(seq_phi, log(density)), collapse = "")
  )
  return (  inla.prior.table )
  }
  
  result = list(seq_phi, density)
  return (result)
  
}



#' The 1d analytic density WCP prior for xi of generalized Pareto distribution
#'
#' @param seq_xi A vector of values of xi parameter. 
#' @param eta User specified parameter of the WCP prior.
#' @param inla_table Should the results be returned as a table to be used in INLA? Default is FALSE.
#'
#' @return A list including densities on seq_xi and a inla table
#' @export
#'
WCP1_1D_GPtail_analytic = function(seq_xi,
                              eta,
                              inla_table = FALSE){
  
  density = numeric()
  for (i in 1:length(seq_xi)){
    xi = seq_xi[i]
    if (xi < 0 | xi >= 1) {
      stop("xi should be a value in [0,1)!")
    }
    density[i] = eta/(1 - xi)^2 * exp( - eta * xi/(1 - xi) )
  }
  
  if (inla_table){  
  # inla table
  inla.prior.table = paste0("table: ",
                            paste(c(seq_xi, density), collapse = "")
  )
  return (inla.prior.table)
  }
  
  result = list(seq_xi, density)
  return (result)
}

#' The 1d analytic density WCP2 prior for tau (1/variance) of Gaussian distribution
#'
#' @param seq_tau A vector of values of tau parameter. 
#' @param eta User specified parameter of the WCP prior.
#' @param inla_table Should the results be returned as a table to be used in INLA? Default is FALSE.
#'
#' @return A list including densities on seq_tau and a inla table
#' @export
#'

WCP2_1D_Gaussian_precision_analytic = function(
                              seq_tau,
                              eta,
                              inla_table = FALSE){

  density = numeric()
  for (i in 1:length(seq_tau)){
    tau = seq_tau[i]
    if (tau < 0 ) {
      stop("tau should be a positive value!")
    }
    density[i] = 0.5 * tau^(-3/2) * eta *exp( -eta * tau^(-1/2))
  }
  if (inla_table){  
  # inla table
  inla.prior.table = paste0("table: ",
                            paste(c(seq_tau, density), collapse = "")
  )
  return (inla.prior.table)
  }
  
  result = list(seq_tau, density )
  return (result)
}

#' The 1d analytic density WCP2 prior for mean of Gaussian distribution
#'
#' @param seq_m A vector of values of mean parameter. 
#' @param eta User specified parameter of the WCP prior.
#' @param inla_table Should the results be returned as a table to be used in INLA? Default is FALSE.
#'
#' @return A list including densities on seq_m and a inla table
#' @export
#'
WCP2_1D_Gaussian_mean_analytic = function(seq_m,
                                              eta,
                                              inla_table = FALSE){
  
  density = numeric()
  for (i in 1:length(seq_m)){
    m = seq_m[i]
    density[i] = 0.5 * eta * exp( - eta * abs(m) )
  }
  if (inla_table){
    # inla table
    inla.prior.table = paste0("table: ",
                              paste(c(seq_m, density), collapse = "")
    )
    return (inla.prior.table)
  }
  
  result = list(seq_m, density )
  return (result)
}

  