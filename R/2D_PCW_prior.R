
## A function to compute the Wasserstein distance between the flexible model and the base model.
## Here we assume that the probability induced by the base model is a Dirac measure.

#' @noRd 

Get_WD = function(flexible_density, 
                  density_arg_name,
                  lower_bound,
                  upper_bound,
                  para_name_flexible, 
                  para_value_flexible,
                  s, 
                  p,
                  args_flexible = list() ){
  
  # flexible_density is the density function of the flexible model
  # density_arg_name is the name of the main parameter of the flexible_density function, for example,
  # in dnorm function, density_arg_name is x.
  # lower_bound and upper_bound specify the range of value of the ramdom variable of flexible model.
  # para_name_flexible is a vector of the names of the parameters (an argument in flexible_density)
  # para_value_flexible is a vector of the values corresponds to the para_name_flexible
  # s is the value where the Dirac measure concentrates at.
  # p is order of the Wasserstein distance
  # args_flexible is a list of additional (fixed) arguments for flexible_density.
  
  
  # update args_flexible
  para_name_flexible = substitute(para_name_flexible)
  i = 1
  for (y in as.list(para_name_flexible)[-1]){
  string  = deparse(y)
  args_flexible[[string]] = para_value_flexible[i]
  i = i + 1
  }
  string = deparse(substitute(density_arg_name)) 
  #print(string)
  flexible_density_unipara = function(x){
      args_flexible[[string]] = x
    return (do.call(flexible_density,args_flexible))
  }
  # check if the density function is a Dirac Delta at s
  # if (flexible_density_unipara(s) == Inf){
  #   # approximate Dirac delta
  #   ddf = function(x){
  #     eps = 1e-6
  #     result = 0.5/sqrt(pi*eps)
  #     result = f*exp(-x^2/(4*eps))
  #     return(result)
  #   }
  # }

  # create integrand
  integrand = function(x){
    temp = flexible_density_unipara(x) * sum((x - s)*(x - s))^(p/2)
    return(abs(temp))
  }
  #integral = integrate(integrand, lower_bound, upper_bound, subdivisions = 200000L,rel.tol = .Machine$double.eps^.05)$value
  integral = cubintegrate(integrand, lower_bound, upper_bound, method = "suave")$integral
  #print(integral)
  return (integral^(1/p))
}
  

################################################# finished testing #############################################



## A function factory to produce a function to compute the Wasserstein distance between the base model and a flexible model.
## The output is a function only require input of parameters of interest.

#' @noRd

WD_function_factory = function(flexible_density, 
                  density_arg_name,
                  lower_bound,
                  upper_bound,
                  para_name_flexible, 
                  s, 
                  p,
                  args_flexible = list() ){
  
  # flexible_density is the density function of the flexible model
  # density_arg_name is the name of the main parameter of the flexible_density function, for example,
  # in dnorm function, density_arg_name is x.
  # lower_bound and upper_bound specify the range of value of the ramdom variable of flexible model.
  # para_name_flexible is a vector of the names of the parameters (an argument in flexible_density)
  # para_value_flexible is a vector of the values corresponds to the para_name_flexible
  # s is the value where the Dirac measure concentrates at.
  # p is order of the Wasserstein distance
  # args_flexible is a list of additional (fixed) arguments for flexible_density.
  
  
  # update args_flexible
  para_name_flexible = substitute(para_name_flexible)
  i = 1
  for (y in as.list(para_name_flexible)[-1]){
  string  = deparse(y)
  args_flexible[[string]] = NA
  i = i + 1
  }
  #print(args_flexible)
  string = deparse(substitute(density_arg_name)) 
  #print(string)
  # the function we should output
  W_2D_func = function(theta1, theta2){
    theta = c(theta1,theta2)
    for (j in 1:length(args_flexible)){
      # update args_flexible
      args_flexible[[j]] = theta[j]
    }
    #print(args_flexible)
    flexible_density_unipara = function(x){
         args_flexible[[string]] = x
         #print(args_flexible)
    return (do.call(flexible_density,args_flexible))
    }
    # create integrand
    integrand = function(x){
    temp = flexible_density_unipara(x) * sum((x - s)*(x - s))^(p/2)
    return(abs(temp))
    }
    # choose a good interval for integration
    a = seq(from = lower_bound,to = upper_bound, length.out = 1000*(upper_bound-lower_bound))
    b = integrand(a)
    #print(tail(b))
    c = which(b>1e-4)
    #print(a[head(c,n = 1)])
    #print(a[tail(c,n = 1)])
    integral = integrate(Vectorize(integrand), a[head(c,n = 1)], a[tail(c,n = 1)], rel.tol=.Machine$double.eps^.05)$value
    #integral = cubintegrate(integrand, lower_bound, upper_bound, method = "suave")$integral
    return (integral^(1/p))
  }
  return(W_2D_func)
  
}


## A function to generate 2d parameter pairs such that the Wasserstein distance among them does not vary much.
## The output is a n by 2 matrix.
## The output will be served as an input to INLA 2D mesh function

#' @noRd


Generating_2D_grid = function(theta_1_init
                              , theta_2_init
                              #,theta_1_range
                              #,theta_2_range
                              , W_function
                              , step_size 
                              , tol.ratio = 0.7
                              #,W_lower
                              #,W_upper
                              ){
    # theta_1_init is a 2d vector, theta_1_init[1] is the lower bound of theta_1 and theta_1_init[2] is the upper bound set by users.
    # The same for theta_2_init

    
    # y records set of values of theta_2 in the resulting grid
    y = c(theta_2_init[1])
    
    # update y
    theta_2_now = theta_2_init[1]
    while (theta_2_now < theta_1_init[2]){
      h = step_size
      W_now = W_func(theta_1_init[1],theta_2_now)
      while ( abs(W_func(theta_1_init[1],theta_2_now+h) - W_now)/W_now > tol.ratio){
        #print(abs(W_func(theta_1_init[1],theta_2_now+h) - W_now)/W_now)
        h = 0.7*h
        #print('h')
        #print(h)
      }
      theta_2_now = theta_2_now + h
      print(theta_2_now)
      if(theta_2_now <= theta_2_init[2]){
      y = c(y, theta_2_now)
      }
    }
    y = c(y, theta_2_init[2])
    # get rid of repeated element
    y = union(y,y)
    print(y)
    

    # coordinates of theta_1 and theta_2
    # theta1_list = c()
    # theta2_list = c()

    # # update x
    # for (i in y){
    #   print('i')
    #   print(i)
    # # append theta_1_init[1]
    # theta_1_now = theta_1_init[1]
    # theta1_list = c(theta1_list, theta_1_now)
    # theta2_list = c(theta2_list, i)
    # while (theta_1_now < theta_1_init[2]){
    #   h = step_size
    #   W_now = W_func(theta_1_now,i)
    #   while ( abs(W_func(theta_1_now+h,i) - W_now) > tol.ratio*W_now){
    #     if (i == y[32]) print(abs(W_func(theta_1_now+h,i) - W_now)/W_now)
    #     #print(abs(W_func(theta_1_now+h,i) - W_now)/W_now)
    #     h = 0.7*h
    #   }
    #   theta_1_now = theta_1_now + h
    #   if(theta_1_now <= theta_1_init[2]){
    #      theta1_list = c(theta1_list, theta_1_now)
    #      theta2_list = c(theta2_list, i)
    #   }
    # }
    # # append theta_1_init[2]
    # theta1_list = c(theta1_list, theta_1_init[2])
    # theta2_list = c(theta2_list, i)
    # }
    # print(length(theta1_list))
    
    # # theta_1 = c()
    # # theta_2 = c()
    # # for (i in x){
    # #   for (j in y){
    # #     theta_1 = c(theta_1,i)
    # #     theta_2 = c(theta_2,j)
    # #   }
    # # }
    # result = as.matrix(rbind(theta1_list, theta2_list))
    # result = t(result)
    # return(result)
}

#' @noRd
grid_2D_generator = function(theta1_low 
                            ,theta1_up
                            ,theta1_step
                            ,theta2_low
                            ,theta2_up
                            ,theta2_step
                            ){
theta1 = seq(from = theta1_low, to = theta1_up, by = theta1_step)
theta2 = seq(from = theta2_low, to = theta2_up, by = theta2_step)
theta1 = c(theta1, theta1_up)
theta1 = union(theta1,theta1)
theta2 = c(theta2, theta2_up)
theta2 = union(theta2,theta2)
grid_coord = expand.grid(theta1,theta2)
grid_coord = as.matrix(grid_coord)
#grid_coord = t(grid_coord)
return(grid_coord)
}

WCP_2D_grid_density = function( theta1_M_low 
                          ,theta1_M_up 
                          ,theta1_M_step 
                          ,theta2_M_low 
                          ,theta2_M_up 
                          ,theta2_M_step 
                          ,W_func
                          ,theta1_D_low 
                          ,theta1_D_up 
                          ,theta1_D_step 
                          ,theta2_D_low 
                          ,theta2_D_up 
                          ,theta2_D_step
                          ,eta
){
print("creating mesh")
# create coordinates for generating fem mesh
domain = grid_2D_generator(theta1_low = theta1_M_low
                            ,theta1_up = theta1_M_up
                            ,theta1_step = theta1_M_step
                            ,theta2_low = theta2_M_low
                            ,theta2_up = theta2_M_up
                            ,theta2_step = theta2_M_step)

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
# make A matrix
A = inla.spde.make.A(mesh, loc = domain )
# weights
print("creating FEM weights")
weights = c()
for (i in 1:dim(mesh$loc)[1]){
  weights = c(weights,W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
}
print("creating coordinates of grid for constructing Jacobian")
# coordinates of grid for construct Jacobian
coord = grid_2D_generator(theta1_low = theta1_D_low 
                            ,theta1_up = theta1_D_up
                            ,theta1_step = theta1_D_step 
                            ,theta2_low = theta2_D_low
                            ,theta2_up = theta2_D_up
                            ,theta2_step = theta2_D_step)

theta1 = seq(from = theta1_D_low, to = theta1_D_up, by = theta1_D_step)
theta1 = c(theta1, theta1_D_up)
theta1 = union(theta1,theta1)
theta2 = seq(from = theta2_D_low, to = theta2_D_up, by = theta2_D_step)
theta2 = c(theta2, theta2_D_up)
theta2 = union(theta2,theta2)
# obtain the Wasserstein distance on each grid point
A_coord = inla.spde.make.A(mesh, loc = coord )
W_value = A_coord %*% weights
# initialize partial arc length on each grid point
parc = numeric(length(W_value))
# initialize total arc length on each grid point
tarc = numeric(length(W_value))
# construct parc and tarc
index = 1:dim(coord)[1]
print("constructing parc and tarc")
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
# initialize the partial derivatives at each grid point
W_partial_theta1 = numeric(dim(coord)[1])
P_partial_theta1 = numeric(dim(coord)[1])
W_partial_theta2 = numeric(dim(coord)[1])
P_partial_theta2 = numeric(dim(coord)[1])
# compute the Jacobian at each coord point
print("computing Jacobian")
for (i in 1:dim(coord)[1]){
    print(i)
    # obtain the parameters
    x = coord[i,1]
    y = coord[i,2]
    W = W_value[i]
    P = parc[i]
    # compute partial derivative of W with respect to theta1
    if (x == theta1_D_low){
        # obtain the right grid point of x
        x_right = x + theta1_D_step
        # obtain the W_value of x_right, y
        temp1 = which(abs( coord[,1] - x_right)<1e-6 )
        temp2 = which(abs( coord[,2] - y)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_right = W_value[temp_index]
        P_right = parc[temp_index]
        # one side differencing
        W_partial_theta1[i] = (W_right - W)/theta1_D_step
        P_partial_theta1[i] = (P_right - P)/theta1_D_step
    } else if (x == theta1_D_up){
        x_left = theta1[length(theta1)-1]
        temp1 = which(abs( coord[,1] - x_left)<1e-6 )
        temp2 = which(abs( coord[,2] - y)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_left = W_value[temp_index]
        P_left = parc[temp_index]
        # one side differencing
        W_partial_theta1[i] = (W - W_left)/theta1_D_step
        P_partial_theta1[i] = (P - P_left)/theta1_D_step
    } else {
        # obtain the right grid point of x
        x_right = x + theta1_D_step
        if(x == theta1[length(theta1)-1]){
            x_right = theta1_D_up
        }
        # obtain the W_value of x_right, y
        temp1 = which(abs( coord[,1] - x_right)<1e-6 )
        temp2 = which(abs( coord[,2] - y)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_right = W_value[temp_index]
        P_right = parc[temp_index]

        # obtain the left grid point of x
        x_left = x - theta1_D_step
        temp1 = which(abs( coord[,1] - x_left)<1e-6 )
        temp2 = which(abs( coord[,2] - y)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_left = W_value[temp_index]
        P_left = parc[temp_index]
        # two side differencing
        W_partial_theta1[i] = (W_right - W_left)/(2 * theta1_D_step)
        P_partial_theta1[i] = (P_right - P_left)/(2 * theta1_D_step)
        }

      # compute partial derivative of W with respect to theta1
    if (y == theta2_D_low){
        # obtain the right grid point of x
        y_up = y + theta2_D_step
        # obtain the W_value of x_right, y
        temp1 = which(abs( coord[,1] - x)<1e-6 )
        temp2 = which(abs( coord[,2] - y_up)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_up = W_value[temp_index]
        P_up = parc[temp_index]
        # one side differencing
        W_partial_theta2[i] = (W_up - W)/theta2_D_step
        P_partial_theta2[i] = (P_up - P)/theta2_D_step
    } else if (y == theta2_D_up){
        y_down = theta2[length(theta2)-1]
        temp1 = which(abs( coord[,1] - x)<1e-6 )
        temp2 = which(abs( coord[,2] - y_down)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_down = W_value[temp_index]
        P_down = parc[temp_index]
        # one side differencing
        W_partial_theta2[i] = (W - W_down)/theta2_D_step
        P_partial_theta2[i] = (P - P_down)/theta2_D_step
    } else {
        # obtain the right grid point of x
        y_up = y + theta2_D_step
        if(y == theta2[length(theta2)-1]){
            y_up = theta2_D_up
        }
        # obtain the W_value of x_right, y
        temp1 = which(abs( coord[,1] - x)<1e-6 )
        temp2 = which(abs( coord[,2] - y_up)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_up = W_value[temp_index]
        P_up = parc[temp_index]

        # obtain the left grid point of x
        y_down = y - theta2_D_step
        temp1 = which(abs( coord[,1] - x)<1e-6 )
        temp2 = which(abs( coord[,2] - y_down)<1e-6 )
        temp_index = intersect(temp1,temp2)
        W_down = W_value[temp_index]
        P_down = parc[temp_index]
        # two side differencing
        W_partial_theta2[i] = (W_up - W_down)/(2 * theta2_D_step)
        P_partial_theta2[i] = (P_up - P_down)/(2 * theta2_D_step)
        }
}
# determinant of density
print("finishing...")
detJ_abs = abs(W_partial_theta1 * P_partial_theta2 - W_partial_theta2 * P_partial_theta1)
approx_WCP_density = eta * detJ_abs * exp(-eta * W_value)/tarc
result = list()
result[[1]] = coord
result[[2]] = approx_WCP_density
result[[3]] = mesh$n
return(result)
}







