library(WCPprior)
library(INLA)
library(excursions)
library(sp)
# using the true Wasserstein distance
W_func = function(sigma, xi){
  true_W = sigma/(1-xi)
  return(true_W)
}

stepsize = c(0.007,0.006,0.005,0.004,0.003,0.002)
L1error = numeric()
for (k in 1:length(stepsize)){
  theta1_low = 0.0001
  theta1_up = 0.3
  theta1_step = stepsize[k]
  theta2_low = 0
  theta2_up = 0.999
  theta2_step = stepsize[k]
  # create coordinates for generating fem mesh
  domain = grid_2D_generator(theta1_low = theta1_low
                             ,theta1_up = theta1_up
                             ,theta1_step = theta1_step
                             ,theta2_low = theta2_low
                             ,theta2_up = theta2_up
                             ,theta2_step = theta2_step)
  # set boundary information for domain as an input to inla.mesh.segment function
  loc.bnd <- matrix(  c(theta1_low, theta2_low, 
                        theta1_up, theta2_low, 
                        theta1_up, theta2_up, 
                        theta1_low,theta2_up)
                      , 4, 2, byrow = TRUE)
  # obtain the segment as an input to inla.mesh.2d
  segm.bnd <- inla.mesh.segment(loc.bnd)
  mesh = inla.mesh.2d(loc = domain
                      ,max.edge = 0.1
                      ,boundary = segm.bnd
  )
  A = inla.spde.make.A( mesh, loc = domain ) 
  weights = c()
  for (i in 1:dim(mesh$loc)[1]){
    #print(i)
    weights = c(weights,W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
  }
  
  theta1_low = 0.001
  theta1_up = 0.3
  theta1_step = 0.003
  theta2_low = 0
  theta2_up = 0.999
  theta2_step = 0.003
  # create coordinates for generating fem mesh
  coord = grid_2D_generator(theta1_low = theta1_low
                            ,theta1_up = theta1_up
                            ,theta1_step = theta1_step
                            ,theta2_low = theta2_low
                            ,theta2_up = theta2_up
                            ,theta2_step = theta2_step)
  A_coord = inla.spde.make.A( mesh, loc = coord) 
  W_value = A_coord%*%weights
  # compute partial and total arc length for the points in the domain
  parc = numeric(length(W_value))
  tarc = numeric(length(W_value))
  W_upper_bound = 0.2
  NA_index = which(W_value > W_upper_bound)
  parc[NA_index] = NA
  tarc[NA_index] = NA
  index = which(W_value <= W_upper_bound)
  while(length(index) > 0){
    #print(index)
    # obtain the current Wasserstein distance
    W = W_value[index[1]]
    print(W)
    print(length(index))
    #print(length(index))
    # obtain all the index of all grid points on this level curve
    grid_level_curve_index = which(W_value == W) 
    
    # obtain the level curve spatial line object with W as the Wasserstein distance
    levelcurve = tricontourmap(mesh, z = weights,
                               tol = 1e-30,
                               levels = c(W))$contour
    # obtain the coordinates of all the discrete points of the level curve
    line_coord = coordinates(levelcurve)[[1]][[1]]
    levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3]
    
    # obtain the partial arc length of the grid point 
    for (i in grid_level_curve_index){
      #print(i)
      # extract the coordinate of the grid point
      grid_point_coordinate = coord[i,] 
      temp1 = which(abs( line_coord[,1] - grid_point_coordinate[1])<1e-6 )
      temp2 = which(abs( line_coord[,2] - grid_point_coordinate[2])<1e-6 )
      temp_index_lc = intersect(temp1,temp2)
      if (length(temp_index_lc) == 0){
        # find the closest x-coord of grid_point_coordinate[1]
        # if (grid_point_coordinate[1] > max(line_coord[,1])){
        
        # }
        coord_index = which(line_coord[,1] == max( line_coord[,1][which(line_coord[,1] < grid_point_coordinate[1])]) )
        x_coord = line_coord[coord_index,1]
        y_coord = line_coord[coord_index,2]
        parc[i] = levelcurve_parc[coord_index] + sqrt((grid_point_coordinate[1] - x_coord)^2 + (grid_point_coordinate[2] - y_coord)^2)
      }else{
        parc[i] = levelcurve_parc[temp_index_lc]
      }
      #print(temp_index_lc)
      tarc[i] = levelcurve_parc[1]
      # remove i from index
      index = index[!index == i]
    }
  }
  theta1 = seq(from = theta1_low, to = theta1_up, by = theta1_step)
  theta1 = c(theta1, theta1_up)
  theta1 = union(theta1,theta1)
  theta2 = seq(from = theta2_low, to = theta2_up, by = theta2_step)
  theta2 = c(theta2, theta2_up)
  theta2 = union(theta2,theta2)
  index = which(W_value <= W_upper_bound)
  W_partial_theta1 = numeric(dim(coord)[1])
  P_partial_theta1 = numeric(dim(coord)[1])
  W_partial_theta1[NA_index] = NA
  P_partial_theta1[NA_index] = NA
  # compute the Jacobian at each coord point
  for (i in index){
    print(i)
    # obtain the parameters
    x = coord[i,1]
    y = coord[i,2]
    W = W_value[i]
    P = parc[i]
    # compute partial derivative of W with respect to theta1
    if (x == theta1_low){
      # obtain the right grid point of x
      x_right = x + theta1_step
      # obtain the W_value of x_right, y
      temp1 = which(abs( coord[,1] - x_right)<1e-6 )
      temp2 = which(abs( coord[,2] - y)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_right = W_value[temp_index]
      P_right = parc[temp_index]
      # one side differencing
      W_partial_theta1[i] = (W_right - W)/theta1_step
      P_partial_theta1[i] = (P_right - P)/theta1_step
    } else if (x == theta1_up){
      x_left = theta1[length(theta1)-1]
      temp1 = which(abs( coord[,1] - x_left)<1e-6 )
      temp2 = which(abs( coord[,2] - y)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_left = W_value[temp_index]
      P_left = parc[temp_index]
      # one side differencing
      W_partial_theta1[i] = (W - W_left)/theta1_step
      P_partial_theta1[i] = (P - P_left)/theta1_step
    } else {
      # obtain the right grid point of x
      x_right = x + theta1_step
      if(x == theta1[length(theta1)-1]){
        x_right = theta1_up
      }
      # obtain the W_value of x_right, y
      temp1 = which(abs( coord[,1] - x_right)<1e-6 )
      temp2 = which(abs( coord[,2] - y)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_right = W_value[temp_index]
      P_right = parc[temp_index]
      
      # obtain the left grid point of x
      x_left = x - theta1_step
      temp1 = which(abs( coord[,1] - x_left)<1e-6 )
      temp2 = which(abs( coord[,2] - y)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_left = W_value[temp_index]
      P_left = parc[temp_index]
      # two side differencing
      W_partial_theta1[i] = (W_right - W_left)/(2 * theta1_step)
      P_partial_theta1[i] = (P_right - P_left)/(2 * theta1_step)
    }
  }
  
  W_partial_theta2 = numeric(dim(coord)[1])
  P_partial_theta2 = numeric(dim(coord)[1])
  W_partial_theta2[NA_index] = NA
  P_partial_theta2[NA_index] = NA
  # compute the Jacobian at each coord point
  for (i in index){
    print(i)
    # obtain the parameters
    x = coord[i,1] 
    y = coord[i,2]
    W = W_value[i]
    P = parc[i]
    # compute partial derivative of W with respect to theta1
    if (y == theta2_low){
      # obtain the right grid point of x
      y_up = y + theta2_step
      # obtain the W_value of x_right, y
      temp1 = which(abs( coord[,1] - x)<1e-6 )
      temp2 = which(abs( coord[,2] - y_up)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_up = W_value[temp_index]
      P_up = parc[temp_index]
      # one side differencing
      W_partial_theta2[i] = (W_up - W)/theta2_step
      P_partial_theta2[i] = (P_up - P)/theta2_step
    } else if (y == theta2_up){
      y_down = theta2[length(theta2)-1]
      temp1 = which(abs( coord[,1] - x)<1e-6 )
      temp2 = which(abs( coord[,2] - y_down)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_down = W_value[temp_index]
      P_down = parc[temp_index]
      # one side differencing
      W_partial_theta2[i] = (W - W_down)/theta2_step
      P_partial_theta2[i] = (P - P_down)/theta2_step
    } else {
      # obtain the right grid point of x
      y_up = y + theta2_step
      if(y == theta2[length(theta2)-1]){
        y_up = theta2_up
      }
      # obtain the W_value of x_right, y
      temp1 = which(abs( coord[,1] - x)<1e-6 )
      temp2 = which(abs( coord[,2] - y_up)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_up = W_value[temp_index]
      P_up = parc[temp_index]
      
      # obtain the left grid point of x
      y_down = y - theta2_step
      temp1 = which(abs( coord[,1] - x)<1e-6 )
      temp2 = which(abs( coord[,2] - y_down)<1e-6 )
      temp_index = intersect(temp1,temp2)
      W_down = W_value[temp_index]
      P_down = parc[temp_index]
      # two side differencing
      W_partial_theta2[i] = (W_up - W_down)/(2 * theta2_step)
      P_partial_theta2[i] = (P_up - P_down)/(2 * theta2_step)
    }
  }
  # compute abs determinant
  detJ_abs = abs(W_partial_theta1 * P_partial_theta2 - W_partial_theta2 * P_partial_theta1)
  eta = 10
  approx_WCP_density = eta * detJ_abs * exp(-eta * W_value)/tarc
  
  W_vector = W_func(coord[,1],coord[,2])
  truedet = (coord[,1]^2+coord[,1]*W_vector^3*(1-coord[,2]))/(W_vector^2*(1-coord[,2])^3*sqrt(W_vector^2+1))
  true_WCP_density = truedet*eta*exp(-eta*W_vector)/sqrt(W_vector^2+1)
  abs_error = abs(approx_WCP_density - true_WCP_density)
  density_index = which((is.na(approx_WCP_density)==F))
  abs_error = abs_error[density_index]
  L1_error = sum(abs_error)*0.1/length(density_index)
  L1error[k] = L1_error
}
#plot
meshsize = c(6336, 8568, 12261, 19076, 33734)#, 75651)
L1_error = c(0.01188243, 0.01100413, 0.01067354, 0.01037834, 0.01009580)#, 0.01036104)
L1_error = 0.5*L1_error
library(ggplot2)
data = cbind(as.vector(meshsize),as.vector(L1_error))
data = data.frame(data)
g1 = ggplot(data, mapping = aes(factor(X1), X2, group = 1)) +
  geom_point(size = 1) + geom_line() +
  labs(x = "Number of mesh nodes", y = expression("Error in L"[1]~"norm") )


g1
ggsave("L1_error_GPareto.png", width = 150,height = 90, unit = 'mm')
