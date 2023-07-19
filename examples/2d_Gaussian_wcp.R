library(INLA)
library(excursions)
library(sp)
# theta1 is mean, theta2 is std
# theta1 = seq(from = -0.15, to = 0.15, by = 0.001)
# theta2 = seq(from = 0.0001, to = 0.15, by = 0.001)
# theta2 = c(theta2, 0.15)
# theta1_list = c()
# theta2_list = c()
# domain = c()
# for (i in theta1){
#   for (j in theta2){
#     theta1_list = c(theta1_list,i)
#     theta2_list = c(theta2_list,j)
#   }
# }
# domain = as.matrix(rbind(theta1_list, theta2_list))
# domain = t(domain)

theta1_low = -0.15
theta1_up = 0.15
theta1_step = 0.005
theta2_low = 0.0001
theta2_up = 0.15
theta2_step = 0.005
# create coordinates for generating fem mesh
domain = grid_2D_generator(theta1_low = theta1_low
                            ,theta1_up = theta1_up
                            ,theta1_step = theta1_step
                            ,theta2_low = theta2_low
                            ,theta2_up = theta2_up
                            ,theta2_step = theta2_step)
# clock-wise can not work
# loc.bnd <- matrix(c(theta1[1], theta2[1]
#                     , theta1[1], tail(theta2,n = 1)
#                     , tail(theta1, n = 1), tail(theta2, n = 1)
#                     , tail(theta1, n = 1), theta2[1]), 4, 2, byrow = TRUE)

# counter clock-wise can work
loc.bnd <- matrix(  c(theta1_low, theta2_low, 
                    theta1_up, theta2_low, 
                    theta1_up, theta2_up, 
                    theta1_low,theta2_up)
                    , 4, 2, byrow = TRUE)
segm.bnd <- inla.mesh.segment(loc.bnd)
mesh = inla.mesh.2d(loc = domain
                    ,max.edge = 0.01
                    ,boundary = segm.bnd
                    )# visualize the mesh
plot(mesh)
points(domain[,1], domain[, 2])
# make A matrix
A = inla.spde.make.A(mesh, loc = domain )
# weights
W_func = WD_function_factory(flexible_density = dnorm, 
                  density_arg_name = x,
                  lower_bound = -100,
                  upper_bound = 100,
                  para_name_flexible = c(mean,sd), 
                  s = 0, 
                  p = 2,
                  args_flexible = list() )
weights = c()
for (i in 1:dim(mesh$loc)[1]){
  weights = c(weights,W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
}
# for (i in 1:nrow(domain)){
#   weights = c(weights,W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
# }

# sup norm error comparing with the true values
true_W = sqrt(domain[,1]^2 + domain[,2]^2)
approx_W = A%*%weights
print(max(abs(approx_W - true_W)))

# construct coord for computing the Jacobian and density
theta1_low = -0.1
theta1_up = 0.1
theta2_low = 0.001
theta2_up = 0.1
theta1_step = 0.003
theta2_step = 0.003
coord = grid_2D_generator(theta1_low = theta1_low 
                            ,theta1_up = theta1_up
                            ,theta1_step = theta1_step 
                            ,theta2_low = theta2_low
                            ,theta2_up = theta2_up
                            ,theta2_step = theta2_step)
theta1 = seq(from = theta1_low, to = theta1_up, by = theta1_step)
theta1 = c(theta1, theta1_up)
theta1 = union(theta1,theta1)
theta2 = seq(from = theta2_low, to = theta2_up, by = theta2_step)
theta2 = c(theta2, theta2_up)
theta2 = union(theta2,theta2)
# theta2 = c(theta2, 0.1)
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
A_coord = inla.spde.make.A(mesh, loc = coord )
W_value = A_coord %*% weights
# sup norm error comparing with true values
true_value = sqrt(coord[,1]^2 + coord[,2]^2)
print(max(abs(true_value - W_value)))



# levelset <- tricontourmap(mesh, z = weights,
#                        levels = c(W_value[1]))$contour
# line_coord = coordinates(levelset)[[1]][[1]]

# compute partial and total arc length for point in coord
parc = numeric(length(W_value))
tarc = numeric(length(W_value))
index = 1:dim(coord)[1]
while(length(index) > 0){
    #print(index)
    # obtain the current Wasserstein distance
    W = W_value[index[1]]
    print(W)
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
        grid_point_coordinate = coord[i,]
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

W_partial_theta1 = numeric(dim(coord)[1])
P_partial_theta1 = numeric(dim(coord)[1])
# compute the Jacobian at each coord point
for (i in 1:dim(coord)[1]){
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
# compute the Jacobian at each coord point
for (i in 1:dim(coord)[1]){
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
eta = 1
approx_WCP_density = eta * detJ_abs * exp(-eta * W_value)/tarc

true_WCP_density = ( 1/sqrt(coord[,1]^2+coord[,2]^2) ) * eta * exp(-eta * sqrt(coord[,1]^2+coord[,2]^2))/pi
abs_error = abs(approx_WCP_density - true_WCP_density)
relative_abs_error = abs_error/true_WCP_density
#l2_log_error = sqrt(sum((log(true_WCP_density) - log(approx_WCP_density))^2))
l2_error = sqrt(sum((true_WCP_density - approx_WCP_density)^2))



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
step_value = c(0.001,0.003,0.005)
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
distance = sqrt(coord[,1]^2+coord[,2]^2)
distance_set[[i]] = 
approx_WCP_density = result[[2]]
true_WCP_density = ( 1/sqrt(coord[,1]^2+coord[,2]^2) ) * eta * exp(-eta * sqrt(coord[,1]^2+coord[,2]^2))/pi
relative_abs_error = abs(approx_WCP_density - true_WCP_density)/abs(true_WCP_density)

max_r_error = c(max_r_error, max(relative_abs_error))
average_r_error = c(average_r_error, sum(relative_abs_error)/length(relative_abs_error))
}

