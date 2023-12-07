#' @noRd 
# function to perform line search along a specified direction

linesearch <- function(theta, s = 1, w, epsilon, W_p, maxit = 1000){
    phi <- atan2(theta[2], theta[1])
    v <- c(cos(phi), sin(phi))
    it <- 0
    while(abs(W_p(theta)-w) > epsilon^2){
        while(W_p(theta) > w){
            theta <- theta - s*v
        }
        while((W_p(theta)<w) && it < maxit){
            theta <- theta + s*v
            it <- it + 1
        }
        s <- s/2
    }
    if(it == maxit){
        stop("Maximum number of iterations reached. Consider increasing 'maxit' to try again. Note that it is possible that W_p does not satisfy the required assumptions.")
    }
    return(theta)
}

#' @noRd 
# function to obtain the approximate level curve

pathfind <- function(phi_l, phi_r, s, epsilon, w, epsilon_tilde, W_p, maxit){
    theta_0 <- linesearch(theta = s*c(cos(phi_l), sin(phi_l)), s = s, w = w, epsilon = epsilon, W_p = W_p, maxit = maxit)
    theta_1 <- linesearch(theta = s*c(cos(phi_r), sin(phi_r)), s = s, w = w, epsilon = epsilon, W_p = W_p, maxit = maxit)
    k <- 1
    theta <- rbind(theta_0,theta_1)
    diff_theta <- diff(theta)
    diff_theta <- sqrt(diff_theta[1]^2 + diff_theta[2]^2)
    while(diff_theta > epsilon_tilde){
        new_theta <- matrix(nrow = 1 + 2^k, ncol = 2)
        for(i in seq(from = 1, to = 2^k + 1, by = 2)){
            new_theta[i,] <- theta[(i-1)/2 + 1,]
        }
        for(i in seq(from = 2, to = 2^k, by = 2)){
            new_theta[i,] <- linesearch(theta = (new_theta[i-1,] + new_theta[i+1,])/2, s = s/2^k, w = w, epsilon = epsilon, W_p = W_p)
        }
        k <- k + 1
        theta <- new_theta
        new_theta <- diff(new_theta)
        diff_theta <- max(apply(new_theta, 1, function(x){sqrt(x[1]^2+x[2]^2)}))
    }
    return(theta)
}

#' region construction
#' @noRd 
# theta_0 -> base model (vertex of the cone)

# W_p <- function(x){sqrt(x[1]^2+x[2]^2)}
# reg <- region_conic(theta_0 = c(0,0), phi_l = 0, phi_r = pi/2, eta = 1, s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p)
# mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)

region_conic <- function(theta_0, phi_l = 0, phi_r, eta, s = 1, epsilon, epsilon_tilde = epsilon, delta, tau = epsilon/1000, W_p, maxit = 1000, return_list = FALSE){
    w_up <- -log(delta/2)/eta
    w_lwr <- -log(1-(delta+tau)/2)/eta + W_p(theta_0 + tau * c(cos((phi_r + phi_l)/2), sin((phi_r+phi_l)/2)))
    theta <- pathfind(phi_l = phi_l, phi_r = phi_r, s = s, epsilon = epsilon, w = w_up, epsilon_tilde =  epsilon_tilde, W_p = W_p, maxit = maxit)
    theta_small <- pathfind(phi_l = phi_l, phi_r = phi_r, s = s, epsilon = epsilon/10, w = w_lwr, epsilon_tilde =  epsilon_tilde/100, W_p = W_p, maxit = maxit)
    theta_small <- theta_small[nrow(theta_small):1,]

    bnd_1 <- theta[1,]
    bnd_2 <- theta[nrow(theta),]

    d_temp <- function(x,y){sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)}

    dist_tmp <- d_temp(bnd_1, theta[2,])
    idx_tmp <- 2
    while(dist_tmp < tau){
        idx_tmp <- idx_tmp + 1
        dist_tmp <- d_temp(bnd_1, theta[idx_tmp,])
    }
    theta <- theta[idx_tmp:nrow(theta),]

    dist_tmp <- d_temp(bnd_2, theta[nrow(theta)-1,])
    idx_tmp <- nrow(theta) - 1
    
    while(dist_tmp < tau){
        idx_tmp <- idx_tmp - 1
        dist_tmp <- d_temp(bnd_2, theta[idx_tmp,])
    }
    theta <- theta[1:idx_tmp,]

    boundary = rbind(theta, theta_small)


    if(return_list){
        return(list(inner = theta_small, outer = theta, boundary = boundary))
    } else{
        return(rbind(theta, theta_small))
    }
}
