library(WCPprior)
library(excursions)
library(ggplot2)
library(foreach)
library(doParallel)
library(viridis)
library(pracma)
library(fmesher)

W_func = function(mean, sd){
  true_W = sqrt(mean^2 + sd^2)
  return(true_W)
}

cutoff = 0.01
mesh_width = 0.5
tau = mesh_width/1000
alpha = 1
eta = 1
W_upper_bound = -log(cutoff/2)/eta
W_lower_bound = -log(1-(cutoff+tau)/2)/eta
# lower and upper bound of theta1 and theta2
L1 = -Inf
U1 = Inf
L2 = 0
U2 = Inf
boundary_path = analytical_boundary_Gaussian(cutoff_out = cutoff/2,
                                             cutoff_inner = cutoff/2 + tau,
                                             eta = eta,
                                             lift = tau,
                                             arc_width = mesh_width^alpha)
boundary_path = matrix(boundary_path, ncol = 2)

result = WCP_2D_Numerical_Density(boundary_path = boundary_path,
                        W_func = W_func,
                        mesh_width = 0.4,
                        alpha = 1,
                        W_lower_bound = W_lower_bound,
                        W_upper_bound = W_upper_bound,
                        eta = 1,
                        L1 = L1, U1 = U1, L2 = L2, U2 = U2,
                        level_curve_type = 'LL',
                        lc_multiplier = 20,
                        NumCores = 3)
density_location = result[[1]]
approx_WCP_density = result[[2]]
mass_lump_C = result[[3]]
A_new = result[[4]]
true_WCP_density = ( 1/sqrt(density_location[,1]^2 + density_location[,2]^2) ) * eta * exp(-eta * sqrt(density_location[,1]^2 + density_location[,2]^2))/pi
TVD_error = 0.5 * (sum(mass_lump_C %*% abs(A_new %*% matrix(approx_WCP_density,ncol = 1) - A_new %*% matrix(true_WCP_density,ncol = 1))) + cutoff)
