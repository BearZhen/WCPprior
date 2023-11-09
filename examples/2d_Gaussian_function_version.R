
library(INLA)
library(excursions)
library(sp)
library(WCPprior)
# function version of the above procedure
W_func = WD_function_factory(flexible_density = dnorm, 
                  density_arg_name = x,
                  lower_bound = -100,
                  upper_bound = 100,
                  para_name_flexible = c(mean,sd), 
                  s = 0, 
                  p = 2,
                  args_flexible = list() )


#max_absr_error = c()
#average_absr_error = c()

error = numeric()
step_value = c(0.005, 0.004, 0.003, 0.002, 0.001)
for (i in 1:length(step_value)){
result = WCP_2D_grid_density( theta1_M_low = -0.15
                          ,theta1_M_up = 0.15
                          ,theta1_M_step = step_value[i]
                          ,theta2_M_low = 0.0001
                          ,theta2_M_up = 0.15
                          ,theta2_M_step = step_value[i]
                          ,W_func = W_func
                          ,theta1_D_low = -0.1
                          ,theta1_D_up = 0.1
                          ,theta1_D_step = 0.003
                          ,theta2_D_low = 0.001
                          ,theta2_D_up = 0.1
                          ,theta2_D_step = 0.003
                          ,eta = 1
)
eta = 1
coord = result[[1]]
#distance = sqrt(coord[,1]^2+coord[,2]^2)
#distance_set[[i]] = distance
approx_WCP_density = result[[2]]
true_WCP_density = ( 1/sqrt(coord[,1]^2+coord[,2]^2) ) * eta * exp(-eta * sqrt(coord[,1]^2+coord[,2]^2))/pi
L1_error = sum(abs(approx_WCP_density - true_WCP_density))*0.2*(0.1-0.003)/dim(coord)[1]
error[i] = L1_error

#max_absr_error = c(max_absr_error, max(relative_abs_error))
#average_absr_error = c(average_absr_error, sum(relative_abs_error)/length(relative_abs_error))
}


library(ggplot2)
data = cbind(as.vector(coord[,1]),as.vector(coord[,2]),as.vector(relative_abs_error))
data = data.frame(data)
g1 = ggplot(data, aes(X1, X2, color = X3))+
geom_point(size = 0.5)+
scale_color_gradient(low="blue", high="red") +
labs(x = expression(theta[1]), y = expression(theta[2])) + 
guides(col = guide_colourbar(title = "RAE"))

g1
ggsave("distribution_of_error.png", width = 150,height = 90, unit = 'mm')





####### test why the error on the right side is big
result = WCP_2D_grid_density( theta1_M_low = -0.15
                          ,theta1_M_up = 0.15
                          ,theta1_M_step = 0.005
                          ,theta2_M_low = 0.0001
                          ,theta2_M_up = 0.15
                          ,theta2_M_step = 0.005
                          ,W_func = W_func
                          ,theta1_D_low = -0.1
                          ,theta1_D_up = 0.11
                          ,theta1_D_step = 0.003
                          ,theta2_D_low = 0.001
                          ,theta2_D_up = 0.1
                          ,theta2_D_step = 0.003
                          ,eta = 1
)
coord = result[[1]]
approx = result[[2]]
true = ( 1/sqrt(coord[,1]^2+coord[,2]^2) ) * eta * exp(-eta * sqrt(coord[,1]^2+coord[,2]^2))/pi
error = abs(approx - true)/abs(true)
data = cbind(as.vector(coord[,1]),as.vector(coord[,2]),as.vector(error))
data = data.frame(data)
g = ggplot(data, aes(X1, X2, color = X3))+
geom_point(size = 0.5)+
scale_color_gradient(low="blue", high="red") +
labs(x = expression(theta[1]), y = expression(theta[2])) + 
guides(col = guide_colourbar(title = "RAE"))
g


theta1_M_low = -0.15
theta1_M_up = 0.15
theta1_M_step = 0.005
theta2_M_low = 0.0001
theta2_M_up = 0.15
theta2_M_step = 0.005

domain = grid_2D_generator(theta1_low = theta1_M_low
                           ,theta1_up = theta1_M_up
                           ,theta1_step = 0.001
                           ,theta2_low = theta2_M_low
                           ,theta2_up = theta2_M_up
                           ,theta2_step = 0.001)
# set up boundary, only counter clock-wise boundary point can work
loc.bnd <- matrix(  c(theta1_M_low, theta2_M_low, 
                      theta1_M_up, theta2_M_low, 
                      theta1_M_up, theta2_M_up, 
                      theta1_M_low,theta2_M_up)
                    , 4, 2, byrow = TRUE)
segm.bnd <- inla.mesh.segment(loc.bnd)
# set up mesh
mesh = inla.mesh.2d(loc = domain
                    ,max.edge = 0.01
                    ,boundary = segm.bnd
)
meshsize = numeric()
meshsize[1] = 1891
meshsize[2] = 2964
meshsize[3] = 5151
meshsize[4] = 11476
meshsize[5] = 45451
L1_error = c(0.002110935, 0.001747103, 0.001734936, 0.001427294, 0.001395999)
L1_error = 0.5*L1_error
### plot number of mesh nodes vs L1 error
library(ggplot2)
data = cbind(as.vector(meshsize),as.vector(L1_error))
data = data.frame(data)
g1 = ggplot(data, mapping = aes(factor(X1), X2, group = 1)) +
  geom_point(size = 1) + geom_line() +
  labs(x = "Number of mesh nodes", y = expression("Error in L"[1]~"norm") )


g1
ggsave("L1_error_Gaussian.png", width = 150,height = 90, unit = 'mm')

