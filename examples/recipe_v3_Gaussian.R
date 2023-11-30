W_func = function(mean, sd){
  true_W = sqrt(mean^2 + sd^2)
  return(true_W)
}

################ construct boundary ########
cutoff = 0.01
mesh_width = 0.005
eta = 46

base_theta1 = 0
base_theta2 = 0
# lower and upper bound of theta1 and theta2
L1 = -Inf
U1 = Inf
L2 = 0
U2 = Inf

# find the upper bound of W based on the cutoff parameter
W_upper_bound = -log(cutoff)/eta

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
    if ( W_func(current_theta1 + mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
      
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
      
    } else if ( W_func(current_theta1, current_theta2 - mesh_width) >= W_upper_bound ){
      current_theta1 = current_theta1 + mesh_width
      current_theta2 = current_theta2 
      # check if the current point is still inside domain, if yes, then update boundary_path
      if (current_theta2 <= U2 & current_theta1 >= L1){
        boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
      }
    } else {
      stop("Fail to construct a boundary path")
    }
  }
  
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
library(fmesher)

#boundary_path = rbind(boundary_path, c(base_theta1, base_theta2))
boundary_path = matrix(boundary_path, ncol = 2)
#mesh = fm_mesh_2d(loc = domain, boundary = fm_segm( boundary_path, is.bnd = TRUE)
#boundary = fm_nonconvex_hull( boundary_path, convex = 0.1)
#    ,max.edge = 0.1
#)
mesh = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width)
plot(mesh,axes = TRUE)


weights = numeric()
for (i in 1:dim(mesh$loc)[1]){
  weights[i] = W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2]))
}

PD_util = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc)
A = PD_util$A
ordering = numeric()
for (i in 1:length(weights)){
ordering[i] = which(abs(A[i,] - 1)<1e-6)
}




############# obtain level curves and filter the non-complete ones
library(excursions)
library(sp)
library(WCPprior)
W_lower_bound = 0.001
parc = numeric()
tarc = numeric()
density_location = numeric()
theta2_low = 0
for (W in seq(from = W_lower_bound, to = W_upper_bound, length = 50)){
  print(W)
  levelcurve = tricontourmap(mesh, z = weights,
                             #tol = 1e-6,
                             levels = W)$contour
  # skip if there is no such level curve
  temp = try(coordinates(levelcurve)[[1]][[1]], silent = FALSE)
  if ('try-error' %in% class(temp)){
    next
  }
  # obtain the coordinates of all the discrete points of the level curve
  line_coord = coordinates(levelcurve)[[1]][[1]]
  # drop the non-complete curve
  if (line_coord[1,2] - theta2_low > 2*mesh_width | line_coord[length(line_coord[,1]),2] - theta2_low > mesh_width){
    next
  }
  # update density location
  density_location = rbind(density_location,line_coord)
  # update 
  levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3]
  parc = c(parc, levelcurve_parc)
  tarc = c(tarc, rep(levelcurve_parc[1], times = length(levelcurve_parc)))
}
#density_location = cbind(density_location,rep(0, times = dim(density_location)[1]))
# construct linear interpolation for parc

PD_util_P = fm_basis(mesh, density_location, derivatives = TRUE)



PD_util_P = fmesher:::fm_evaluator_mesh_2d(mesh, density_location, derivatives = TRUE)
P_PD_x = PD_util_P$dx %*% matrix(parc, ncol = 1)
P_PD_x = as.vector(P_PD_x)
# compute partial derivatives of W_func with respect to y, sigma
P_PD_y = PD_util_P$dy %*% matrix(parc, ncol = 1)
P_PD_y = as.vector(P_PD_y)







