# try 2d interpolation for 2d Generalized Pareto W distance
theta1_low = 0.001
theta1_up = 0.3
theta1_step = 0.005
theta2_low = 0
theta2_up = 0.999
theta2_step = 0.005
theta1 = seq(from = theta1_low, to = theta1_up, by = theta1_step)
theta1 = c(theta1, theta1_up)
theta1 = union(theta1,theta1)
theta2 = seq(from = theta2_low, to = theta2_up, by = theta2_step)
theta2 = c(theta2, theta2_up)
theta2 = union(theta2,theta2)
grid <- meshgrid(theta1, theta2)
W = grid$X/(1 - grid$Y)

theta1_low = 0.001
theta1_up = 0.25
theta2_low = 0
theta2_up = 0.999
theta1_step = 0.006
theta2_step = 0.006
coord = grid_2D_generator(theta1_low = theta1_low 
                            ,theta1_up = theta1_up
                            ,theta1_step = theta1_step 
                            ,theta2_low = theta2_low
                            ,theta2_up = theta2_up
                            ,theta2_step = theta2_step)
# theta1p = seq(from = theta1_low, to = theta1_up, by = theta1_step)
# theta1p = c(theta1p, theta1_up)
# theta1p = union(theta1p,theta1p)
# theta2p = seq(from = theta2_low, to = theta2_up, by = theta2_step)
# theta2p = c(theta2p, theta2_up)
# theta2p = union(theta2p,theta2p)

method <- "linear"
Wp <- interp2(theta1, theta2, W, coord[,1], coord[,2], method)
W_func = function(sigma, xi){
  true_W = sigma/(1-xi)
  return(true_W)
}

Wt = W_func(coord[,1],coord[,2])
abs_error = abs(Wt-Wp)

method <- "linear"
Wp <- interp2(theta1, theta2, W, coord[,1], coord[,2], method)
abs_error = abs(Wt-Wp)

# now try with barycentric Lagrange interpolation
theta1_low = 0.001
theta1_up = 0.25
theta2_low = 0
theta2_up = 0.999
theta1_step = 0.006
theta2_step = 0.006
theta1p = seq(from = theta1_low, to = theta1_up, by = theta1_step)
theta1p = c(theta1p, theta1_up)
theta1p = union(theta1p,theta1p)
theta2p = seq(from = theta2_low, to = theta2_up, by = theta2_step)
theta2p = c(theta2p, theta2_up)
theta2p = union(theta2p,theta2p)

Wp = barylag2d(t(W), theta1, theta2, theta1p, theta2p)
Wgrid = meshgrid(theta1p,theta2p)
Wt = t(Wgrid$X/(1 - Wgrid$Y))
abs_error = abs(Wt-Wp)
