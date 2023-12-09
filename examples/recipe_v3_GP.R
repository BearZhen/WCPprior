library(excursions)
library(sp)
library(WCPprior)
library(akima)
library(ggplot2)
library(fmesher)


W_func = function(sigma, xi){
  true_W = sigma/(1-xi)
  return(true_W)
}

################ construct boundary ########
cutoff = 0.01
mesh_width = 0.001
eta = 46

base_theta1 = 0
base_theta2 = 0
# lower and upper bound of theta1 and theta2
L1 = 0
U1 = Inf
L2 = 0
U2 = 1

# find the upper bound of W based on the cutoff parameter
W_upper_bound = -log(cutoff/2)/eta
W_lower_bound = -log(1-cutoff/2)/eta
Z_star = NA
boundary_path = c()


# theta2 is bounded, theta1 is one side bounded
# we assume that base_theta1 = L1, base_theta2 = L2
if (U2 - L2 < Inf & L1 > -Inf & U1 == Inf){
  # find a direction to search
  N = L1
  current_W = W_func(N, (L2+U2)/2 )
  
  # what if W_func(L1, (L2+U2)/2 ) >= W_upper_bound? Then we should directly construct boundary path
  
  # if not, we do binary search to find Z_star
  while ( current_W < W_upper_bound ){
    N = N + 1
    current_W = W_func(N, (L2+U2)/2 )
  }
  BS_interval_upper = N
  BS_interval_lower = N - 1
  BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
  one_step_back_W = W_func( BS_interval_upper - mesh_width, (L2+U2)/2 )
  while ( one_step_back_W >= W_upper_bound ){
    if ( W_func( BS_interval_mid, (L2+U2)/2 ) < W_upper_bound){
      BS_interval_lower = BS_interval_mid
      BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
    } else{
      BS_interval_upper = BS_interval_mid
      BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
      # update
      one_step_back_W = W_func( BS_interval_upper - mesh_width, (L2+U2)/2 )
      #current_W = W_func( BS_interval_upper, (L2+U2)/2 )
    }
  }
  Z_star = BS_interval_upper
  
  
  
  
  # cosntruct the boundary path when theta2 towards L2
  starting_theta1 = Z_star
  starting_theta2 = (L2+U2)/2 
  current_theta1 = starting_theta1
  current_theta2 = starting_theta2
  boundary_path = rbind( boundary_path, c(starting_theta1, starting_theta2) )
  
  while ( current_theta2 > L2 ){
    # find the next point
    if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
      current_theta1 = current_theta1 
      current_theta2 = current_theta2 - mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 >= L2){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1 + mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
      
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 - mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 >= L2){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1 + mesh_width, current_theta2 ) >= W_upper_bound ){
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 
      boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      
    } else if ( W_func(current_theta1 + mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 + mesh_width
      boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      
    } else {
      print(c(current_theta1, current_theta2))
      stop("Fail to construct a boundary path")
      
    }
  }
  boundary_path[,1] = rev(boundary_path[,1])
  boundary_path[,2] = rev(boundary_path[,2])
  # cosntruct the boundary path when theta2 towards U2
  current_theta1 = starting_theta1
  current_theta2 = starting_theta2
  
  while ( current_theta2 < U2 ){
    # find the next point
    if ( W_func(current_theta1 - mesh_width, current_theta2 + mesh_width) >= W_upper_bound ){
      current_theta1 = current_theta1 - mesh_width
      current_theta2 = current_theta2 + mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 <= U2 & current_theta1 >= L1){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1, current_theta2 + mesh_width) >= W_upper_bound ){
      current_theta1 = current_theta1  
      current_theta2 = current_theta2 + mesh_width 
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 <= U2 & current_theta1 >= L1){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
    } else {
      print(c(current_theta1, current_theta2))
      stop("Fail to construct a boundary path")
      
    }
  }
  boundary_path = rbind( boundary_path, c( min(boundary_path[,1]), min(boundary_path[,2]) ) )
  
} else if (L1 == -Inf & U1 == Inf & L2 > -Inf & L2 < Inf & U2 == Inf) {
  # theta2 is one side bounded, theta1 is two side unbounded
  # we assume that base_theta2 = L2
  
  # find Z_star
  N = L2
  current_W = W_func(base_theta1, N)
  while ( current_W < W_upper_bound ){
    N = N + 1
    current_W = W_func(base_theta1, N )
  }
  BS_interval_upper = N
  BS_interval_lower = N - 1
  BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
  one_step_back_W = W_func( base_theta1, BS_interval_upper - mesh_width )
  while ( one_step_back_W >= W_upper_bound ){
    if ( W_func( base_theta1, BS_interval_mid ) < W_upper_bound){
      BS_interval_lower = BS_interval_mid
      BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
    } else{
      BS_interval_upper = BS_interval_mid
      BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
      # update
      one_step_back_W = W_func( base_theta1, BS_interval_upper - mesh_width )
      #current_W = W_func( BS_interval_upper, (L2+U2)/2 )
    }
  }
  Z_star = BS_interval_upper
  # cosntruct the boundary path when theta1 towards U1
  starting_theta2 = Z_star
  starting_theta1 = base_theta1 
  current_theta2 = starting_theta2
  current_theta1 = starting_theta1
  boundary_path = rbind( boundary_path, c(starting_theta1, starting_theta2) )
  while ( current_theta2 > L2 ){
    # find the next point
    if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
      current_theta1 = current_theta1 
      current_theta2 = current_theta2 - mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 >= L2){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1 + mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
      
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 - mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 >= L2){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1 + mesh_width, current_theta2 ) >= W_upper_bound ){
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 
      boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      
    } else if ( W_func(current_theta1 + mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 + mesh_width
      boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      
    } else {
      stop("Fail to construct a boundary path")
    }
  }
  
  boundary_path[,1] = rev(boundary_path[,1])
  boundary_path[,2] = rev(boundary_path[,2])
  
  current_theta2 = starting_theta2
  current_theta1 = starting_theta1
  while ( current_theta2 > L2 ){
    # find the next point
    if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
      current_theta1 = current_theta1 
      current_theta2 = current_theta2 - mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 >= L2){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1 - mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
      
      current_theta1 = current_theta1 - mesh_width
      current_theta2 = current_theta2 - mesh_width
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 >= L2){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
      
    } else if ( W_func(current_theta1 - mesh_width, current_theta2 ) >= W_upper_bound ){
      current_theta1 = current_theta1 - mesh_width
      current_theta2 = current_theta2 
      boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      
    } else if ( W_func(current_theta1 - mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
      current_theta1 = current_theta1 - mesh_width
      current_theta2 = current_theta2 + mesh_width
      boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      
    } else {
      stop("Fail to construct a boundary path")
    }
  }
  
}


################## construct mesh

#boundary_path = rbind(boundary_path, c(base_theta1, base_theta2))
boundary_path = matrix(boundary_path, ncol = 2)
#mesh = fm_mesh_2d(loc = domain, boundary = fm_segm( boundary_path, is.bnd = TRUE)
#boundary = fm_nonconvex_hull( boundary_path, convex = 0.1)
#    ,max.edge = 0.1
#)
mesh = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width)
plot(mesh,axes = TRUE)



