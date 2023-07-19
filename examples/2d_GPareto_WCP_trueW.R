# density function for the generalized Pareto distribution
# GPdensity = function(x,sigma, xi){
#     result = sigma^(-1)*(1+xi*x/sigma)^(-1-1/xi)
#     return(result)
# }

#x = seq(from = 0, to = 200,by = 0.01)
#plot(x,GPdensity(x, 0.0001, 0.9999))

library(WCPprior)
library(INLA)
# using the true Wasserstein distance
W_func = function(sigma, xi){
  true_W = sigma/(1-xi)
  return(true_W)
}

# sigma = 0.0001
# xi = 0.9999
# true_W = sigma/(1-xi)
# print(true_W)
# print(W_func(sigma,xi))

W_end = 0.2
# theta1 is sigma, theta2 is xi
theta1_low = 0.001
theta1_up = 0.3
theta1_step = 0.01
theta2_low = 0
theta2_up = 0.999
theta2_step = 0.01
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
                    )# visualize the mesh
plot(mesh)
points(domain[,1], domain[, 2])
# make A matrix
A = inla.spde.make.A( mesh, loc = domain ) 
weights = c()
for (i in 1:dim(mesh$loc)[1]){
  print(i)
  weights = c(weights,W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
}
# for (i in 1:nrow(domain)){
#   weights = c(weights,W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
# }

# sup norm error comparing with the true values
true_W = W_func(domain[,1], domain[,2])
approx_W = A%*%weights
W_abs_error = abs(approx_W - true_W)
print(max(abs(approx_W - true_W)))





