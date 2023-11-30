#rm(list = ls())


#library(INLA)
library(fmesher)
library(ggplot2)
library(excursions)
library(sp)
library(WCPprior)
# create a inla mesh

theta1_low = -0.15
theta1_up = 0.15
theta1_step = 0.004
theta2_low = 0.0001
theta2_up = 0.15
theta2_step = 0.004
# create coordinates for generating fem mesh
domain = grid_2D_generator(theta1_low = theta1_low
                             ,theta1_up = theta1_up
                             ,theta1_step = theta1_step
                             ,theta2_low = theta2_low
                             ,theta2_up = theta2_up
                             ,theta2_step = theta2_step)

# loc.bnd <- matrix(  c(theta1_low, theta2_low, 
#                         theta1_up, theta2_low, 
#                         theta1_up, theta2_up, 
#                         theta1_low,theta2_up)
#                       , 4, 2, byrow = TRUE)

  # obtain the segment as an input to inla.mesh.2d
#segm.bnd <- inla.mesh.segment(loc.bnd)
# mesh1 = inla.mesh.2d(loc = domain
#                       ,max.edge = 0.1
#                       ,boundary = segm.bnd
# )
# plot(mesh1)

# Now, create the same mesh with fmesher 
mesh = fm_mesh_2d(
    loc = domain
    ,boundary = fm_segm(rbind(c(theta1_low, theta2_low), 
                            c(theta1_up, theta2_low), 
                            c(theta1_up, theta2_up), 
                            c(theta1_low,theta2_up)), is.bnd = TRUE)
    ,max.edge = 0.1
)
plot(mesh, axes = TRUE)


# create a Wasserstein function
W_func = WD_function_factory(flexible_density = dnorm, 
                             density_arg_name = x,
                             lower_bound = -100,
                             upper_bound = 100,
                             para_name_flexible = c(mean,sd), 
                             s = 0, 
                             p = 2,
                             args_flexible = list() )

# weights1 = numeric()
# for (i in 1:dim(mesh1$loc)[1]){
#   weights1[i] = W_func(as.numeric(mesh1$loc[i,1]), as.numeric(mesh1$loc[i,2]))
# }

weights = numeric()
for (i in 1:dim(mesh$loc)[1]){
  weights[i] = W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2]))
}

true_W = sqrt(mesh$loc[,1]^2 + mesh$loc[,2]^2)
#proj1 = fm_evaluator(mesh1, loc = mesh1$loc)
#approx_W1 = fm_evaluate(proj1, weights1)
#proj2 = fm_evaluator(mesh2, loc = mesh2$loc)
#PD_util = fm_basis(mesh2, loc = mesh2$loc, weights = weights2, derivatives = TRUE)
PD_util = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc)
#approx_W2 = fm_evaluate(proj2, weights2)
approx_W = PD_util$A %*% matrix(weights, ncol = 1)
approx_W  = as.vector(approx_W)

#print(max(abs(approx_W1 - true_W)))
#print(mean(abs(approx_W - true_W)))
#print(mean(abs(weights - true_W)))
#print(mean(abs(weights - approx_W)))
#print(max(abs(weights - approx_W)))
# obtain an object that can help us to compute partial derivatives of 
#PD_util = fm_basis(mesh2, loc = mesh2$loc, weights = weights2, derivatives = TRUE)
# compute partial derivatives of W_func with respect to x, mean
#W_PD_x = PD_util_W$dx %*% matrix(weights,ncol = 1)
#W_PD_x = as.vector(W_PD_x)
# compute partial derivatives of W_func with respect to y, sigma
#W_PD_y = PD_util_W$dy %*% matrix(weights,ncol = 1)
#W_PD_y = as.vector(W_PD_y)
# true dW/dx
#true_dWdx = mesh$loc[,1]/sqrt(mesh$loc[,1]^2 + mesh$loc[,2]^2)
#true_dWdy = mesh$loc[,2]/sqrt(mesh$loc[,1]^2 + mesh$loc[,2]^2)
#print(mean(abs(W_PD_x - true_dWdx)))
#print(mean(abs(W_PD_y - true_dWdy)))

#coord = mesh$loc
parc = numeric(length(approx_W))
tarc = numeric(length(approx_W))
W_upper_bound = 0.1
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

#index = 1:length(approx_W)
while(length(index) > 0){
  #print(index)
  # obtain the current Wasserstein distance
  W = approx_W[index[1]]
  print(W)
  #print(length(index))
  # obtain all the index of all grid points on this level curve
  grid_level_curve_index = which(approx_W == W) 
  # obtain the level curve spatial line object with W as the Wasserstein distance
  levelcurve = tricontourmap(mesh, z = approx_W,
                             #tol = 1e-6,
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
  #levelcurve_parc = numeric(dim(line_coord)[1])
  #cum_partial_arc_length = 0
  # obtain the partial arc length for all the discrete points on the level curve
  # for (i in 2:dim(line_coord)[1]){
  #   levelcurve_parc[dim(line_coord)[1] - i + 1] = cum_partial_arc_length + sqrt( 
  #     (line_coord[dim(line_coord)[1]-i+2,1]-line_coord[dim(line_coord)[1]-i+1,1])^2 
  #     + (line_coord[dim(line_coord)[1]-i+2,2]-line_coord[dim(line_coord)[1]-i+1,2])^2 )
  #   cum_partial_arc_length = levelcurve_parc[dim(line_coord)[1] - i + 1]
  # }
  # obtain the partial arc length of the grid point 
  for (i in grid_level_curve_index){
    print(i)
    # extract the coordinate of the grid point
    grid_point_coordinate = mesh$loc[i,]
    temp1 = which(abs( line_coord[,1] - grid_point_coordinate[1])<1e-6 )
    temp2 = which(abs( line_coord[,2] - grid_point_coordinate[2])<1e-6 )
    temp_index_lc = intersect(temp1,temp2)
    if (length(temp_index_lc) == 0){
      # find the closest x-coord of grid_point_coordinate[1]
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
true_dPdx = mesh$loc[index,1]*(pi - acos(mesh$loc[index,1]/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2)))/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2) + mesh$loc[index,2]/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2)
true_dPdy = mesh$loc[index,2]*(pi - acos(mesh$loc[index,1]/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2)))/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2) - mesh$loc[index,1]/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2)
print(mean(abs(P_PD_x - true_dPdx), na.rm = TRUE))
print(mean(abs(P_PD_y - true_dPdy), na.rm = TRUE))

PD_util_W = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc[index,], derivatives = TRUE)
W_PD_x = PD_util_W$dx %*% matrix(weights,ncol = 1)
W_PD_x = as.vector(W_PD_x)
# compute partial derivatives of W_func with respect to y, sigma
W_PD_y = PD_util_W$dy %*% matrix(weights,ncol = 1)
W_PD_y = as.vector(W_PD_y)
true_dWdx = mesh$loc[index,1]/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2)
true_dWdy = mesh$loc[index,2]/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2)
print(mean(abs(W_PD_x - true_dWdx)))
print(mean(abs(W_PD_y - true_dWdy)))

approx_detJ_abs = abs(W_PD_x * P_PD_y - W_PD_y * P_PD_x)
true_detJ_abs = abs(true_dWdx*true_dPdy - true_dWdy*true_dPdx)
abs_error_det = mean(abs(approx_detJ_abs - true_detJ_abs), na.rm = TRUE)
eta = 1
approx_WCP_density = eta * approx_detJ_abs * exp(-eta * approx_W[index])/tarc[index]
true_WCP_density = ( 1/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2) ) * eta * exp(-eta * sqrt(mesh$loc[index,1]^2+mesh$loc[index,2]^2))/pi
abs_error = abs(approx_WCP_density - true_WCP_density)
L1_error = sum(abs_error, na.rm = TRUE)*pi*(0.1)^2/sum(!is.na(abs_error))









data = cbind(as.vector(mesh$loc[index,1]),as.vector(mesh$loc[index,2]),as.vector(approx_WCP_density))
data = data.frame(data)
g = ggplot(data, aes(X1, X2, color = X3))+
  geom_point(size = 0.5)+
  scale_color_gradient(low="blue", high="red") +
  labs(x = expression(theta[1]), y = expression(theta[2])) + 
  guides(col = guide_colourbar(title = "Density"))+
  theme(legend.position = "top",
        legend.key.height = unit(0.1, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=5),#change legend text font size
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) 
g
ggsave("approximated_density_Gaussian.png", width = 100,height = 60, unit = 'mm')

gm = ggplot() +
  geom_fm(data = mesh) +
  labs(x = expression(theta[1]), y = expression(theta[2]))+
  theme(
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) 

gm
ggsave("mesh.png", width = 100,height = 60, unit = 'mm')

library(ggpubr)
figure = ggarrange(gm, g, # remove axis labels from plots
                   ncol = 2, nrow = 1,
                   legend = "top",
                   align = "hv"
)
figure
ggsave(filename="mesh_and_density.png", width = 120, height = 60, unit = 'mm')

# plot
TVD = c(0.0023815282, 0.0016114837, 0.0016063054, 0.0011707505, 0.0006571597)
stepsize = c(0.005,0.004,0.003,0.002,0.001)
library(ggplot2)
data = cbind(as.vector(stepsize),as.vector(TVD))
data = data.frame(data)
g1 = ggplot(data, mapping = aes(X1, X2)) +
  geom_point(size = 1) + geom_line() +
  labs(x = "mesh width", y = expression("Total variation distance") ) 




g1
#ggsave("TVD_error_Gaussian.png", width = 100,height = 60, unit = 'mm')











