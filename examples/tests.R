
library(cubature)



##################################### test Get_WD with normal distribution #####################################
result = Get_WD(flexible_density = dnorm, 
                  density_arg_name = x,
                  lower_bound = -Inf,
                  upper_bound = Inf,
                  para_name_flexible = c(mean,sd), 
                  para_value_flexible = c(1,1),
                  s = 0, 
                  p = 2,
                  args_flexible = list() )
print(result)




######################################## test WD_function_factory ###########################################
## generate a function to compute Wasserstein distance between a Dirac measure and Gaussian with mean and std.
W_func = WD_function_factory(flexible_density = dnorm, 
                  density_arg_name = x,
                  lower_bound = -Inf,
                  upper_bound = Inf,
                  para_name_flexible = c(mean,sd), 
                  s = 0, 
                  p = 2,
                  args_flexible = list() )
W_func(0,0)
W_func(1,1)
W_func(-0.1,0.001)
theta_1 = seq(from = -0.1, to = 0.1, by = 0.01)
theta_2 = seq(from = 0.01, to = 0.1, by = 0.01)
for (i in theta_1){
  for (j in theta_2){
    print(c(i,j))
    print(W_func(i,j))
    
  }
}
################################################# finished testing ##########################################



# ## a function to compute two dimensional PCW prior density
# PCW_2D_density = function(theta1_range,
#                           theta2_range,
#                           theta1_step,
#                           theta2_step,
#                           W_function ){

# # construct
# theta1_grid = seq(theta1_range[1],theta1_range[2], by = theta1_step)
# theta1_grid = append(theta1_grid, theta1_range[2])
# theta1_grid = union(theta1_grid, theta1_grid)

# theta2_grid = seq(theta2_range[1], theta2_range[2], by = theta2_step)
# theta2_grid = append(theta2_grid, theta2_range[2])
# theta2_grid = union(theta2_grid, theta2_grid)

# W_matrix = outer(theta1_grid, theta2_grid, W_function)

# return(W_matrix)
# }



# W_matrix = PCW_2D_density(theta1_range = c(-0.1,0.1),
#                           theta2_range = c(0,0.1),
#                           theta1_step = 0.01,
#                           theta2_step = 0.01,
#                           Wasserstein_function = W_func)



# coord = Generating_2D_grid(theta_1_init = c(-0.1,0.1)
#                           , theta_2_init = c(0.001,0.01)
#                           , W_function = W_func
#                           , step_size = 0.01
#                           , tol.ratio = 0.7
#                           )

######################### to be continued from above ##################


# library(INLA)
# theta1 = seq(from = -0.1, to = 0.1, by = 0.01)
# theta2 = seq(from = 0.001, to = 0.1, by = 0.01)

# theta1 = seq(from = 0, to = 1, by =0.1)
# theta2 = seq(from = 0, to = 1, by =0.1)
# theta1_list = c()
# theta2_list = c()
# coord = c()
# for (i in theta1){
#   for (j in theta2){
#     theta1_list = c(theta1_list,i)
#     theta2_list = c(theta2_list,j)
#   }
# }
# coord = as.matrix(rbind(theta1_list, theta2_list))
# coord = t(coord)
# # finite element mesh

# loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
# segm.bnd <- inla.mesh.segment(loc.bnd)

# mesh = inla.mesh.2d(loc = coord
#                     ,max.edge = 0.1,
#                     boundary = segm.bnd
#                     #,offset = c(0,0.1)
#                     )
# # plot mesh and mesh nodes
# plot(mesh)
# points(coord[,1], coord[, 2])
# # make A matrix
# A = inla.spde.make.A(mesh, loc = coord )
# # test
# #A = inla.spde.make.A(mesh, loc = matrix(c(-0.1,0.011),ncol = 2, byrow = T) )
# # make the vector of finite element weights
# weights = c()
# for (i in 1:nrow(coord)){
#   weights = c(weights,sqrt(as.numeric(coord[i,1])^2 + as.numeric(coord[i,2])^2))
# }
# weight_vector = numeric(ncol(A))
# weight_vector[colSums(A)>0] = weights
# new_coord = matrix(c())







