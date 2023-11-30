W_func = function(mean, sd){
  true_W = sqrt(mean^2 + sd^2)
  return(true_W)
}


cutoff = 0.01
mesh_width = 0.003
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

library(ggplot2)
point_path = c()
Z_star = NA
N = L2
current_W = W_func(base_theta1, N)
point_path = rbind(point_path, c(base_theta1, N))
while ( current_W < W_upper_bound ){
  N = N + 1
  current_W = W_func(base_theta1, N )
  point_path = rbind(point_path, c(base_theta1, N))
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
    point_path = rbind(point_path, c(base_theta1, BS_interval_upper))
  }
}
Z_star = BS_interval_upper
point_path = data.frame(point_path)
g = ggplot(point_path, aes(x= X1, y= X2)) + geom_point(size = 1.5) 
                                                                
for (i in 1:(length(point_path$X1) - 1)){
  g = g + geom_curve(x = point_path[i,1], y = point_path[i,2],
                     xend = point_path[i+1,1], yend = point_path[i+1,2],
                     color = 2,
                     curvature = -0.2,
                     arrow = arrow(length = unit(2, "mm")))
}


# cosntruct the boundary path when theta1 towards U1
right_path = c()

starting_theta2 = Z_star
starting_theta1 = base_theta1 

current_theta2 = starting_theta2
current_theta1 = starting_theta1
right_path = rbind( right_path, c(starting_theta1, starting_theta2) )
while ( current_theta2 > L2 ){
  # find the next point
  if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
    current_theta1 = current_theta1 
    current_theta2 = current_theta2 - mesh_width
    # check if the current point is still inside domain, if yes, then update boundary_path
    if (current_theta2 >= L2){
      right_path = rbind( right_path, c(current_theta1, current_theta2) )
    }
    
  } else if ( W_func(current_theta1 + mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
    
    current_theta1 = current_theta1 + mesh_width
    current_theta2 = current_theta2 - mesh_width
    # check if the current point is still inside domain, if yes, then update boundary_path
    if (current_theta2 >= L2){
      right_path = rbind( right_path, c(current_theta1, current_theta2) )
    }
    
  } else if ( W_func(current_theta1 + mesh_width, current_theta2 ) >= W_upper_bound ){
    current_theta1 = current_theta1 + mesh_width
    current_theta2 = current_theta2 
    right_path = rbind( right_path, c(current_theta1, current_theta2) )
    
  } else if ( W_func(current_theta1 + mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
    current_theta1 = current_theta1 + mesh_width
    current_theta2 = current_theta2 + mesh_width
    right_path = rbind( right_path, c(current_theta1, current_theta2) )
    
  } else {
    stop("Fail to construct a boundary path")
  }
}

right_path = data.frame(right_path)
g_boundary = ggplot() + geom_point(right_path, mapping = aes(x= X1, y= X2), size = 0.5) 
g_boundary
for (i in 1:(length(right_path$X1) - 1)){
  g_boundary = g_boundary + geom_curve(right_path, mapping = aes(x= X1, y= X2), x = right_path[i,1], y = right_path[i,2],
                     xend = right_path[i+1,1], yend = right_path[i+1,2],
                     color = 2,
                     curvature = -0.5,
                     arrow = arrow(length = unit(1, "mm")))
}
g_boundary



left_path = c()
current_theta2 = starting_theta2
current_theta1 = starting_theta1
left_path = rbind(left_path, c(current_theta1, current_theta2))
while ( current_theta2 > L2 ){
  # find the next point
  if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
    current_theta1 = current_theta1 
    current_theta2 = current_theta2 - mesh_width
    
    # check if the current point is still inside domain, if yes, then update boundary_path
    if (current_theta2 >= L2){
      left_path = rbind( left_path, c(current_theta1, current_theta2) )
    }
    
  } else if ( W_func(current_theta1 - mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
    
    current_theta1 = current_theta1 - mesh_width
    current_theta2 = current_theta2 - mesh_width
    # check if the current point is still inside domain, if yes, then update boundary_path
    if (current_theta2 >= L2){
      left_path = rbind( left_path, c(current_theta1, current_theta2) )
    }
    
  } else if ( W_func(current_theta1 - mesh_width, current_theta2 ) >= W_upper_bound ){
    current_theta1 = current_theta1 - mesh_width
    current_theta2 = current_theta2 
    left_path = rbind( left_path, c(current_theta1, current_theta2) )
    
  } else if ( W_func(current_theta1 - mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
    current_theta1 = current_theta1 - mesh_width
    current_theta2 = current_theta2 + mesh_width
    left_path = rbind( left_path, c(current_theta1, current_theta2) )
    
  } else {
    stop("Fail to construct a boundary path")
  }
}

left_path = data.frame(left_path)
g_boundary = g_boundary + geom_point(left_path,mapping = aes(x= X1, y= X2), size = 0.5) 
g_boundary
for (i in 1:(length(right_path$X1) - 1)){
  g_boundary = g_boundary + geom_curve(left_path, mapping = aes(x= X1, y= X2), x = left_path[i,1], y = left_path[i,2],
                                       xend = left_path[i+1,1], yend = left_path[i+1,2],
                                       color = 2,
                                       curvature = 0.5,
                                       arrow = arrow(length = unit(1, "mm")))
}
g_boundary

theta1_low = min(boundary_path[,1])
theta1_up = max(boundary_path[,1])
theta1_step = mesh_width
theta2_low = min(boundary_path[,2])
theta2_up = max(boundary_path[,2])
theta2_step = mesh_width
# create coordinates for generating fem mesh
domain = grid_2D_generator(theta1_low = theta1_low
                           ,theta1_up = theta1_up
                           ,theta1_step = theta1_step
                           ,theta2_low = theta2_low
                           ,theta2_up = theta2_up
                           ,theta2_step = theta2_step)

############### construct finite element mesh within the boundary path
library(fmesher)
#boundary_path = rbind(boundary_path, c(base_theta1, base_theta2))
boundary_path = matrix(boundary_path, ncol = 2)
mesh = fm_mesh_2d(loc = domain, boundary = fm_segm( boundary_path, is.bnd = TRUE)
                  #boundary = fm_nonconvex_hull( boundary_path, convex = 0.1)
                  ,max.edge = 0.1
)
plot(mesh,axes = TRUE)