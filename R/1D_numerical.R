
# function to compute gradient
my.grad = function(func, x){
  eps = 1e-6
  return ((func(x + eps)-func(x))/eps)
  
  #require(numDeriv)
  #return (grad(func, x))
}

#' The 1d numerical density WCP prior 
#'
#' @param base_theta Base model theta.
#' @param L Lower bound of theta.
#' @param L_included If L should be included in the domain, a TRUE or FALSE value
#' @param U Upper bound of theta.
#' @param U_included If U should be included in the domain, a TRUE or FALSE value
#' @param W_func A function that returns the Wasserstein distance given theta.
#' @param eta User specified parameter of the WCP prior.
#' @param cutoff1 Cutoff parameter that decide upper bound of the Wasserstein distance on the right hand side of base_theta.
#' @param cutoff2 Cutoff parameter that decide upper bound of the Wasserstein distance on the left hand side of base_theta.
#' @param mesh_width Mesh width to be used. If not used, mesh_n will be used.
#' @param mesh_n Number of mesh nodes to be used. Will not be used if mesh_width is given.
#' @param inla_table Should the results be returned as a table to be used in INLA? Default is FALSE.
#' @return A list of density locations and densities evaluated on those locations.
#' @export
#'
WCP_1D_Numerical_Density= function(W_func,
                                   eta,
                                   base_theta,
                                   L, 
                                   U,
                                   mesh_width = NULL,
                                   mesh_n = 100,
                                   cutoff1 = 0.01,
                                   cutoff2 = NULL,
                                   L_included = FALSE,
                                   U_included = FALSE,
                                   inla_table = FALSE){
  
   
  # when the domain of theta is one-sided 
  if (L == -Inf & base_theta == U){
    if (base_theta == Inf | base_theta == -Inf) {
      stop("The base model parameter should be a finite value!")
    }
    if (U == Inf | U == -Inf) {
      stop("The upper bound of theta should be a finite value!")
    }
    # determine an upper bound for W
    W_upper_bound = -log(cutoff1)/eta
    ## binary search
    N = U
    current_W = W_func(N)
    # if not, we do binary search to find Z_star
    while ( current_W < W_upper_bound ){
      N = N - 1
      current_W = W_func(N)
    }
    BS_interval_lower = N
    BS_interval_upper = N + 1
    BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
    one_step_back_W = W_func( BS_interval_lower + mesh_width )
    while ( one_step_back_W >= W_upper_bound ){
      if ( W_func( BS_interval_mid ) < W_upper_bound){
        BS_interval_upper = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
      } else{
        BS_interval_lower = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
        # update
        one_step_back_W = W_func( BS_interval_lower + mesh_width )
      }
    }
    Z_star = BS_interval_lower
    # create mesh
    mesh = seq(from = Z_star,to = U, by = mesh_width)
    if (U_included == TRUE){
      mesh = union(mesh, U)
    } else {
      mesh = mesh[!mesh == U]
    }
    
    
  } else if (U == Inf & base_theta == L) { # another case of one-sided domain
    if (base_theta == Inf | base_theta == -Inf) {
      stop("The base model parameter should be a finite value!")
    }
    if (L == Inf | L == -Inf) {
      stop("The lower bound of theta should be a finite value!")
    }
    # determine an upper bound for W
    W_upper_bound = -log(cutoff1)/eta
    ## binary search
    N = L
    current_W = W_func(N)
    # if not, we do binary search to find Z_star
    while ( current_W < W_upper_bound ){
      N = N + 1
      current_W = W_func(N)
    }
    BS_interval_upper = N
    BS_interval_lower = N - 1
    BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
    one_step_back_W = W_func( BS_interval_upper - mesh_width )
    while ( one_step_back_W >= W_upper_bound ){
      if ( W_func( BS_interval_mid ) < W_upper_bound){
        BS_interval_lower = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
      } else{
        BS_interval_upper = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
        # update
        one_step_back_W = W_func( BS_interval_upper - mesh_width )
        
      }
    }
    Z_star = BS_interval_upper
    if (L_included == TRUE){
      mesh = seq(from = L,to = Z_star, by = mesh_width)
    } else{
      mesh = seq(from = L + mesh_width, to = Z_star, by = mesh_width )
    }
    mesh = union(mesh, Z_star)
    
  } else if (L == -Inf & U == Inf) { # When the domain is double sided
    if (base_theta == Inf | base_theta == -Inf) {
      stop("The base model parameter should be a finite value!")
    }
    if (is.null(cutoff2) == 1) {
      stop("User should provide cutoff2!")
    }
    W_upper_bound_right = -log(cutoff1)/eta
    W_upper_bound_left = -log(cutoff2)/eta
    ## binary search for lower bound of the interval
    N = base_theta
    current_W = W_func(N)
    # if not, we do binary search to find Z_star
    while ( current_W < W_upper_bound ){
      N = N - 1
      current_W = W_func(N)
    }
    BS_interval_lower = N
    BS_interval_upper = N + 1
    BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
    one_step_back_W = W_func( BS_interval_lower + mesh_width )
    while ( one_step_back_W >= W_upper_bound ){
      if ( W_func( BS_interval_mid ) < W_upper_bound){
        BS_interval_upper = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
      } else{
        BS_interval_lower = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
        # update
        one_step_back_W = W_func( BS_interval_lower + mesh_width )
      }
    }
    Z_star_lower = BS_interval_lower
    
    ## binary search for upper bound of the interval
    N = base_theta
    current_W = W_func(N)
    # if not, we do binary search to find Z_star
    while ( current_W < W_upper_bound ){
      N = N + 1
      current_W = W_func(N)
    }
    BS_interval_upper = N
    BS_interval_lower = N - 1
    BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
    one_step_back_W = W_func( BS_interval_upper - mesh_width )
    while ( one_step_back_W >= W_upper_bound ){
      if ( W_func( BS_interval_mid ) < W_upper_bound){
        BS_interval_lower = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
      } else{
        BS_interval_upper = BS_interval_mid
        BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
        # update
        one_step_back_W = W_func( BS_interval_upper - mesh_width )
        
      }
    }
    Z_star_upper = BS_interval_upper
    # create mesh
    mesh = seq(from = Z_star_lower,to = Z_star_upper, by = mesh_width)
    # add Z_star to the mesh if it is not included
    mesh = union(mesh, Z_star_upper)
    
    
  } else { # when the domain is bounded
    if (base_theta == Inf | base_theta == -Inf) {
      stop("The base model parameter should be a finite value!")
    }
    if (L == Inf | L == -Inf) {
      stop("The lower bound of theta should be a finite value!")
    }
    if (U == Inf | U == -Inf) {
      stop("The upper bound of theta should be a finite value!")
    }
    
    if (L_included == TRUE){
      mesh = seq(from = L,to = U, by = mesh_width)
    } else{
      mesh = seq(from = L + mesh_width, to = U, by = mesh_width )
    }
    
    if (U_included == TRUE){
      mesh = union(mesh, U)
    } else {
      mesh = mesh[!mesh == U]
    }
    
  } 
  # obtain the Wasserstein-2 distance W(xi)
  W = c()
  location = c()
  absJ = c()
  for (i in mesh){
    temp = try(W_func(i), silent = FALSE)
    if ('try-error' %in% class(temp)){
      print(i)
      next
    }
    temp = try(my.grad(W_func, i), silent = FALSE)
    if ('try-error' %in% class(temp)){
      print(i)
      next
    }
    
    location = c(location, i)
    W = c(W,W_func(i))
    J = my.grad(W_func, i)
    absJ = c(absJ, abs(J))
  }
  # obtain dw/dxi
  
  # return the value of density
  density_discrete = absJ*eta*exp(-eta*W)
  # for INLA

  if(inla_table){
    prior_table = paste0("table: ",
                         paste(c(location, density_discrete), collapse = " "))
    return(prior_table)
  }
  
  return(cbind(location, density_discrete))
}





