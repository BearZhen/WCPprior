my.grad = function(func, x){
  eps = 1e-6
  return ((func(x + eps)-func(x))/eps)
  
  #require(numDeriv)
  #return (grad(func, x))
}


# 1d WCP prior density
wcp_density_1d = function(base_theta, L, U, W_func, eta, cutoff = NULL, mesh_width){
  # determine an upper bound for W

  

  if (L == -Inf & base_theta == U){
    
  } else if (U == Inf & base_theta == L) {
    W_upper_bound = -log(cutoff)/eta
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
    # create mesh
    mesh = seq(from = L,to = Z_star, by = mesh_width)
    # add Z_star to the mesh if it is not included
    mesh = union(mesh, Z_star)
    #density_discrete = my.grad(W_fun,mesh)*eta*exp(-eta*W_func(mesh))
    
    
  } else if (U < Inf & base_theta == L){
    # create mesh
    mesh = seq(from = L,to = U, by = mesh_width)
    # add Z_star to the mesh if it is not included
    #mesh = union(mesh, U)
  } else {
    stop("Fail to fina an upper bound for the Wasserstein distance.")
  }
  
  
  
  
  # obtain the Wasserstein-2 distance W(xi)
  W = W_func(mesh)
  # obtain dw/dxi
  dwdxi = my.grad(W_func, mesh)
  # return the value of density
  density_discrete = abs(dwdxi)*eta*exp(-eta*W)
  return (density_discrete)
}


W_func = function(xi){
  return(xi/(1-xi))
}

density = wcp_density_1d(base_theta = 0, L = 0, U = 1, W_func = W_func, eta = 4.6, cutoff = NULL, mesh_width = 0.003)



