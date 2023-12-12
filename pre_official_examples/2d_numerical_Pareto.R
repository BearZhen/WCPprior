rm(list = ls())
library(WCPprior)
library(excursions)
library(ggplot2)
library(foreach)
library(doParallel)
library(viridis)
library(pracma)
library(fmesher)

W_GP = function(x){
  true_W = x[1]/(1-x[2])
  return(true_W)
}

cutoff = 0.01
#mesh_width = 0.0005
alpha = 1
eta = 20
# lower and upper bound of theta1 and theta2
# L1 = 0
# U1 = Inf
# L2 = 0
# U2 = 1
# W_upper_bound = -log(cutoff/2)/eta
# boundary_path = analytical_boundary_GP(cutoff = cutoff,
#                                        eta = eta,
#                                        bound_away_from_one = 0.001)
# boundary_path = matrix(boundary_path, ncol = 2)
# 
# result = WCP_2D_Numerical_Density(boundary_path = boundary_path,
#                         W_func = W_GP,
#                         mesh_width = 0.05,
#                         alpha = 1,
#                         W_lower_bound = 0,
#                         W_upper_bound = W_upper_bound,
#                         eta = 20,
#                         L2 = L2, U2 = U2,
#                         level_curve_type = 'LU',
#                         lc_multiplier = 20,
#                         NumCores = 2)

result = WCP_2D_Numerical_Density(W_GP,
                                  eta = 20,
                                  mesh_width = 0.04, 
                                  alpha = 1,
                                  tau = 0.001,
                                  cutoff = 0.01, 
                                  region = list(type = 'strip', lower_boundary = 0, upper_boundary = 1, direction = 'positive', base_theta = c(0,0) ),
                                  lc_multiplier = 20,
                                  parallel = TRUE,
                                  NumCores = 3)
#density_location = result[[1]]
#approx_WCP_density = result[[2]]
#mass_lump_C = result[[3]]
#A_new = result[[4]]
#W_distance = result[[5]]
true_WCP_density = eta/(1 - result[[1]][,2])* exp(-eta * result[[5]])
TVD_error = 0.5 * (sum(result[[3]] %*% abs(result[[4]] %*% matrix(result[[2]],ncol = 1) - result[[4]] %*% matrix(true_WCP_density,ncol = 1))) + cutoff)
