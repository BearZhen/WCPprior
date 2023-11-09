library(WCPprior)
library(fmesher)
#library(INLA)
library(excursions)
library(sp)
# using the true Wasserstein distance
W_func = function(sigma, xi){
  true_W = sigma/(1-xi)
  return(true_W)
}

stepsize = c(0.006,0.005,0.004,0.003,0.002)
TVD = numeric()
abs_det_error = numeric()
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
  mesh = fm_mesh_2d(
    loc = domain
    ,boundary = fm_segm(rbind(c(theta1_low, theta2_low), 
                              c(theta1_up, theta2_low), 
                              c(theta1_up, theta2_up), 
                              c(theta1_low,theta2_up)), is.bnd = TRUE)
    ,max.edge = 0.1
  )
  plot(mesh, axes = TRUE)
  weights = c()
  for (i in 1:dim(mesh$loc)[1]){
    #print(i)
    weights[i] = W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2]))
  }
  
  PD_util = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc)
  #approx_W2 = fm_evaluate(proj2, weights2)
  approx_W = PD_util$A %*% matrix(weights, ncol = 1)
  approx_W  = as.vector(approx_W)
  
  parc = numeric(length(approx_W))
  tarc = numeric(length(approx_W))
  W_upper_bound = 0.2
  W_lower_bound = 0.001
  NA_index = which(approx_W > W_upper_bound)
  parc[NA_index] = NA
  tarc[NA_index] = NA
  NA_index = which(approx_W < W_lower_bound)
  parc[NA_index] = NA
  tarc[NA_index] = NA
  index1 = which(approx_W <= W_upper_bound)
  index2 = which(approx_W >= W_lower_bound)
  index = intersect(index1, index2)
  
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
    
    temp = try(coordinates(levelcurve)[[1]][[1]], silent = FALSE)
    if ('try-error' %in% class(temp)){
      for (i in grid_level_curve_index){
        parc[i] = NA
        tarc[i] = NA
        index = index[!index == i]
      }
      next
    }
    
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
  
  # construct finite element interpolation for parc
  index1 = which(approx_W <= W_upper_bound)
  index2 = which(approx_W >= W_lower_bound)
  index = intersect(index1, index2)
  PD_util_P = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc[index,], derivatives = TRUE)
  P_PD_x = PD_util_P$dx %*% matrix(parc, ncol = 1)
  P_PD_x = as.vector(P_PD_x)
  # compute partial derivatives of W_func with respect to y, sigma
  P_PD_y = PD_util_P$dy %*% matrix(parc, ncol = 1)
  P_PD_y = as.vector(P_PD_y)
  
  PD_util_W = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc[index,], derivatives = TRUE)
  W_PD_x = PD_util_W$dx %*% matrix(weights,ncol = 1)
  W_PD_x = as.vector(W_PD_x)
  # compute partial derivatives of W_func with respect to y, sigma
  W_PD_y = PD_util_W$dy %*% matrix(weights,ncol = 1)
  W_PD_y = as.vector(W_PD_y)
  
  # compute abs determinant
  approx_detJ_abs = abs(W_PD_x * P_PD_y - W_PD_y * P_PD_x)
  eta = 10
  approx_WCP_density = eta * approx_detJ_abs * exp(-eta * approx_W[index])/tarc[index]
  
  W_vector = W_func(mesh$loc[index,1], mesh$loc[index,2])
  truedet = (mesh$loc[index,1]^2 + mesh$loc[index,1]*W_vector^3*(1-mesh$loc[index,2]))/(W_vector^2*(1-mesh$loc[index,2])^3*sqrt(W_vector^2+1))
  abs_det_error[k] = mean(abs(approx_detJ_abs - truedet),na.rm = TRUE)
  true_WCP_density = truedet*eta*exp(-eta*W_vector)/sqrt(W_vector^2+1)
  abs_error = abs(approx_WCP_density - true_WCP_density)
  L1_error = sum(abs_error, na.rm = TRUE)*0.1/sum(!is.na(abs_error))
  TVD[k] = 0.5*L1_error
}



TVD = c(0.04181517, 0.03891213, 0.03580960, 0.03304342, 0.02569051)
stepsize = c(0.006,0.005,0.004,0.003,0.002)
library(ggplot2)
data = cbind(as.vector(stepsize),as.vector(TVD))
data = data.frame(data)
g1 = ggplot(data, mapping = aes(X1, X2)) +
  geom_point(size = 1) + geom_line() +
  labs(x = "Mesh width", y = expression("Total variation distance") ) 



g1
ggsave("TVD_error_GPareto.png", width = 100,height = 60, unit = 'mm')


