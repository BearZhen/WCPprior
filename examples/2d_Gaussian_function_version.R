
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


max_absr_error = c()
average_absr_error = c()
distance_set = list()
error_set = list()
step_value = c(0.005,0.003,0.001)
for (i in 2:length(step_value)){
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
distance = sqrt(coord[,1]^2+coord[,2]^2)
distance_set[[i]] = distance
approx_WCP_density = result[[2]]
true_WCP_density = ( 1/sqrt(coord[,1]^2+coord[,2]^2) ) * eta * exp(-eta * sqrt(coord[,1]^2+coord[,2]^2))/pi
relative_abs_error = abs(approx_WCP_density - true_WCP_density)/abs(true_WCP_density)
error_set[[i]] = relative_abs_error
max_absr_error = c(max_absr_error, max(relative_abs_error))
average_absr_error = c(average_absr_error, sum(relative_abs_error)/length(relative_abs_error))
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
