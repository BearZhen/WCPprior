
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
    #integral = integrate(Vectorize(integrand), lower_bound, upper_bound, rel.tol=1e-15)$value
    integral = cubintegrate(integrand, lower_bound, upper_bound, method = "suave")$integral
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


