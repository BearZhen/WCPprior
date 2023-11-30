

# cutoff, probability mass cutoff parameter of the (truncated) exponential distribution of the Wasserstein distance.
# mesh_width
# eta, tail mass control parameter
# base_theta1 and base_theta2 are base parameter values for theta1 and theta2


W_Gaussian = function(mean, sd){
  true_W = sqrt(mean^2 + sd^2)
  return(true_W)
}

cutoff = 0.01
mesh_width = 0.003
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

## find an upper bound for theta2 when theta1 = base_thata1 based on W_upper_bound by binary search.

N_star = NA
max.interation = 1000000
parameter_index = 0

# if base_theta2 = L2
if (base_theta2 == L2){
    counter = 0
    search_step = mesh_width
    # search N_star of theta2
    while (counter <= max.interation){
     
     current_theta2 = base_theta2 + search_step
     # break the loop if current_theta2 is outside the domain
     if (current_theta2 > U2){
        break
     }
     # break the loop if we find N_star
     current_W = W_Gaussian(base_theta1, current_theta2)
     if (current_W >= W_upper_bound){
        N_star = current_theta2
        parameter_index = 2
        break
     }
     # update search_step
     search_step = search_step + mesh_width
     # update counter
     counter = counter + 1
    }

    # search N_star of theta1 if theta2 does not work
    if ( is.na(N_star) ){
        # in the case that base_theta1 == L1
       if(base_theta1 == L1){
         search_step = mesh_width
         while (counter <= max.interation){
     
            current_theta1 = base_theta1 + search_step
            # break the loop if current_theta1 is outside the domain
            if (current_theta1 > U2){
             break
            }
            # break the loop if we find N_star
            current_W = W_Gaussian(current_theta1, base_theta2)
            if (current_W >= W_upper_bound){
             N_star = current_theta1
                parameter_index = 1
                break
            }
            # update search_step
            search_step = search_step + mesh_width
            # update counter
            counter = counter + 1
            }
       }
       

    }



}

# if both search in theta1 and theta2 does not work
if ( is.na(N_star) ){
 stop("Fail to find a good (theta1, base_theta2) or (base_theta1, theta2) to meet W_upper_bound .")
}

## construct boundary
# if we found N_star of theta2
if (parameter_index == 2){
    
    # first case, base_theta1 = L1
    if (base_theta1 == L1){
       
    }


}





