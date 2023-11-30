# cutoff, probability mass cutoff parameter of the (truncated) exponential distribution of the Wasserstein distance.
# mesh_width
# eta, tail mass control parameter
# base_theta1 and base_theta2 are base parameter values for theta1 and theta2

library(fmesher)
library(ggplot2)
library(excursions)
library(sp)
library(WCPprior)

W_func = function(mean, sd){
  true_W = sqrt(mean^2 + sd^2)
  return(true_W)
}

cutoff = 0.01
mesh_width = 0.005
eta = 46

base_theta1 = 0
base_theta2 = 0
# lower and upper bound of theta1 and theta2
L1 = -Inf
U1 = Inf
L2 = 0
U2 = Inf

# find the upper bound of W based on the cutoff parameter
W_upper_bound = -log(cutoff)/eta

Z_star = NA
boundary_path = c()
# theta2 is bounded, theta1 is one side bounded
# we assume that base_theta1 = L1, base_theta2 = L2
if (U2 - L2 < Inf & L1 > -Inf & U1 == Inf){
   # find a direction to search
   N = L1
   current_W = W_func(N, (L2+U2)/2 )
   
   # what if W_func(L1, (L2+U2)/2 ) >= W_upper_bound? Then we should directly construct boundary path
   
   # if not, we do binary search to find Z_star
   while ( current_W < W_upper_bound ){
    N = N + 1
    current_W = W_func(N, (L2+U2)/2 )
   }
   BS_interval_upper = N
   BS_interval_lower = N - 1
   BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
   one_step_back_W = W_func( BS_interval_upper - mesh_width, (L2+U2)/2 )
   while ( one_step_back_W >= W_upper_bound ){
         if ( W_func( BS_interval_mid, (L2+U2)/2 ) < W_upper_bound){
            BS_interval_lower = BS_interval_mid
            BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
         } else{
            BS_interval_upper = BS_interval_mid
            BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
            # update
            one_step_back_W = W_func( BS_interval_upper - mesh_width, (L2+U2)/2 )
            #current_W = W_func( BS_interval_upper, (L2+U2)/2 )
         }
   }
   Z_star = BS_interval_upper
   
   
   
   
   # cosntruct the boundary path when theta2 towards L2
   starting_theta1 = Z_star
   starting_theta2 = (L2+U2)/2 
   current_theta1 = starting_theta1
   current_theta2 = starting_theta2
   boundary_path = rbind( boundary_path, c(starting_theta1, starting_theta2) )

   while ( current_theta2 > L2 ){
         # find the next point
         if ( W_func(current_theta1 + mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
           
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 - mesh_width
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 >= L2){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }
            
         } else if ( W_func(current_theta1 + mesh_width, current_theta2 ) >= W_upper_bound ){
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 
            boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )

         } else if ( W_func(current_theta1 + mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 + mesh_width
            boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )

         } else {
            stop("Fail to construct a boundary path")
         }
   }

   # cosntruct the boundary path when theta2 towards U2
   current_theta1 = starting_theta1
   current_theta2 = starting_theta2
   while ( current_theta2 < U2 ){
         # find the next point
         if ( W_func(current_theta1 - mesh_width, current_theta2 + mesh_width) >= W_upper_bound ){
            current_theta1 = current_theta1 - mesh_width
            current_theta2 = current_theta2 + mesh_width
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 <= U2 & current_theta1 >= L1){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }
            
        } else if ( W_func(current_theta1, current_theta2 - mesh_width) >= W_upper_bound ){
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 <= U2 & current_theta1 >= L1){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }
        } else {
            stop("Fail to construct a boundary path")
        }
   }

} else if (L1 == -Inf & U1 == Inf & L2 > -Inf & L2 < Inf & U2 == Inf) {
# theta2 is one side bounded, theta1 is two side unbounded
# we assume that base_theta2 = L2

# find Z_star
    N = L2
    current_W = W_func(base_theta1, N)
    while ( current_W < W_upper_bound ){
    N = N + 1
    current_W = W_func(base_theta1, N )
   }
   BS_interval_upper = N
   BS_interval_lower = N - 1
   BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
   one_step_back_W = W_func( base_theta1, BS_interval_upper - mesh_width )
   while ( one_step_back_W >= W_upper_bound ){
         if ( W_func( base_theta1, BS_interval_mid ) < W_upper_bound){
            BS_interval_lower = BS_interval_mid
            BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
         } else{
            BS_interval_upper = BS_interval_mid
            BS_interval_mid = (BS_interval_upper + BS_interval_lower)/2
            # update
            one_step_back_W = W_func( base_theta1, BS_interval_upper - mesh_width )
            #current_W = W_func( BS_interval_upper, (L2+U2)/2 )
         }
   }
   Z_star = BS_interval_upper
   # cosntruct the boundary path when theta1 towards U1
   starting_theta2 = Z_star
   starting_theta1 = base_theta1 
   current_theta2 = starting_theta2
   current_theta1 = starting_theta1
   boundary_path = rbind( boundary_path, c(starting_theta1, starting_theta2) )
   while ( current_theta2 > L2 ){
         # find the next point
         if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
            current_theta1 = current_theta1 
            current_theta2 = current_theta2 - mesh_width
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 >= L2){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }

         } else if ( W_func(current_theta1 + mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
           
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 - mesh_width
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 >= L2){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }
            
         } else if ( W_func(current_theta1 + mesh_width, current_theta2 ) >= W_upper_bound ){
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 
            boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )

         } else if ( W_func(current_theta1 + mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
            current_theta1 = current_theta1 + mesh_width
            current_theta2 = current_theta2 + mesh_width
            boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )

         } else {
            stop("Fail to construct a boundary path")
         }
   }

   boundary_path[,1] = rev(boundary_path[,1])
   boundary_path[,2] = rev(boundary_path[,2])

   current_theta2 = starting_theta2
   current_theta1 = starting_theta1
   while ( current_theta2 > L2 ){
         # find the next point
         if ( W_func(current_theta1 , current_theta2 - mesh_width) >= W_upper_bound ) {
            current_theta1 = current_theta1 
            current_theta2 = current_theta2 - mesh_width
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 >= L2){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }

         } else if ( W_func(current_theta1 - mesh_width, current_theta2 - mesh_width) >= W_upper_bound ){
           
            current_theta1 = current_theta1 - mesh_width
            current_theta2 = current_theta2 - mesh_width
            # check if the current point is still inside domain, if yes, then update boundary_path
            if (current_theta2 >= L2){
               boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )
            }
            
         } else if ( W_func(current_theta1 - mesh_width, current_theta2 ) >= W_upper_bound ){
            current_theta1 = current_theta1 - mesh_width
            current_theta2 = current_theta2 
            boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )

         } else if ( W_func(current_theta1 - mesh_width, current_theta2 + mesh_width ) >= W_upper_bound  ){
            current_theta1 = current_theta1 - mesh_width
            current_theta2 = current_theta2 + mesh_width
            boundary_path = rbind( boundary_path, c(current_theta1, current_theta2) )

         } else {
            stop("Fail to construct a boundary path")
         }
   }




}



plot(boundary_path[,1], boundary_path[,2])

theta1_low = min(boundary_path[,1])
theta1_up = max(boundary_path[,1])
theta1_step = mesh_width
theta2_low = min(boundary_path[,2])
theta2_up = max(boundary_path[,2])
theta2_step = mesh_width
# create coordinates for generating fem mesh
#domain = grid_2D_generator(theta1_low = theta1_low
#                             ,theta1_up = theta1_up
#                             ,theta1_step = theta1_step
#                             ,theta2_low = theta2_low
#                             ,theta2_up = theta2_up
#                             ,theta2_step = theta2_step)

############### construct finite element mesh within the boundary path
library(fmesher)
#boundary_path = rbind(boundary_path, c(base_theta1, base_theta2))
boundary_path = matrix(boundary_path, ncol = 2)
#mesh = fm_mesh_2d(loc = domain, boundary = fm_segm( boundary_path, is.bnd = TRUE)
     #boundary = fm_nonconvex_hull( boundary_path, convex = 0.1)
#    ,max.edge = 0.1
#)
mesh = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = 0.004)
plot(mesh,axes = TRUE)

weights = numeric()
for (i in 1:dim(mesh$loc)[1]){
  weights[i] = W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2]))
}

PD_util = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc)
#idx = mesh$idx$loc
#idx = idx[!is.na(idx)]
#idx = unique(idx)
A = PD_util$A



approx_W = A %*% matrix(weights, ncol = 1)
approx_W  = as.vector(approx_W)

#ordering = numeric()
#for (i in 1:length(weights)){
  #ordering[i] = which(abs(A[i,] - 1)<1e-6)
#}

#reordered_weights = weights[ordering]

#parc = numeric(length(weights))
parc = numeric(length(approx_W))
tarc = numeric(length(approx_W))
#tarc = numeric(length(weights))
NA_index = which(approx_W > W_upper_bound)
#NA_index = which(weights > W_upper_bound)
parc[NA_index] = NA
tarc[NA_index] = NA
index = which(approx_W <= W_upper_bound)
#index = which(weights <= W_lower_bound)
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
  # skip if there is nothing in levelcurve
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
  # abandon incomplete curves, this part need to be adjusted according to the geometry of the computation domain
  if (line_coord[1,2] - theta2_low > 0.005 | line_coord[length(line_coord[,1]),2] - theta2_low > 0.005){
    for (i in grid_level_curve_index){
      parc[i] = NA
      tarc[i] = NA
      index = index[!index == i]
    } 
    next
  }
  
  
  levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3]
  # obtain the partial arc length of the grid point 
  for (i in grid_level_curve_index){
    print(i)
    # extract the coordinate of the grid point
    #grid_point_coordinate = mesh$loc[ordering[i],]
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


# construct linear interpolation for parc
index = which(!is.na(parc))
PD_util_P = fmesher:::fm_evaluator_mesh_2d(mesh, mesh$loc[index,], derivatives = TRUE)
P_PD_x = PD_util_P$dx %*% matrix(parc, ncol = 1)
P_PD_x = as.vector(P_PD_x)
index_na = which(is.na(parc))
#parc[index_na] = NA
#P_PD_x_test = PD_util_P$dx %*% matrix(parc, ncol = 1)
#P_PD_x_test = as.vector(P_PD_x_test)
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
print(abs_error_det)

abs_exponential_error = mean(abs(exp(-eta * approx_W[index])/tarc[index] - exp(-eta * sqrt(mesh$loc[index,1]^2+mesh$loc[index,2]^2))/pi   ), na.rm = TRUE )


approx_WCP_density = eta * approx_detJ_abs * exp(-eta * approx_W[index])/tarc[index]
true_WCP_density = ( 1/sqrt(mesh$loc[index,1]^2 + mesh$loc[index,2]^2) ) * eta * exp(-eta * sqrt(mesh$loc[index,1]^2+mesh$loc[index,2]^2))/pi
abs_error = abs(approx_WCP_density - true_WCP_density)

#true_WCP_density_1 = true_detJ_abs* eta * exp(-eta * sqrt(mesh$loc[index,1]^2+mesh$loc[index,2]^2))/pi
#abs_error_1 = abs(approx_WCP_density - true_WCP_density_1)



TVD_error = 0.5 * sum(abs_error, na.rm = TRUE)*pi*(0.1)^2/sum(!is.na(abs_error))

TVD_error2 = 0.5 * sum(abs_error, na.rm = TRUE)*mesh_width^2

#TVD_error_1 = 0.5 * sum(abs_error_1, na.rm = TRUE)*pi*(0.1)^2/sum(!is.na(abs_error_1))



data = cbind(as.vector(mesh$loc[index,1]),as.vector(mesh$loc[index,2]),as.vector(abs_error))
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

