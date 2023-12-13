#' @noRd 
# function to perform line search along a specified direction


#' @noRd 

cut_boundary <- function(theta, tau){
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
    return(theta)
}

linesearch_conic <- function(theta, s = 1, w, epsilon, W_p, maxit = 1000){
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
        it <- 0
        s <- s/2
    }
    if(it == maxit){
        stop("Maximum number of iterations reached. Consider increasing 'maxit' to try again. Note that it is possible that W_p does not satisfy the required assumptions.")
    }
    return(theta)
}


#' @noRd 
# function to obtain the approximate level curve

pathfind_conic <- function(phi_l, phi_r, s, epsilon, w, epsilon_tilde, W_p, maxit, buff = 0){
    theta_0 <- linesearch_conic(theta = s*c(cos(phi_l + buff), sin(phi_l + buff)), s = s, w = w, epsilon = epsilon, W_p = W_p, maxit = maxit)
    theta_1 <- linesearch_conic(theta = s*c(cos(phi_r - buff), sin(phi_r - buff)), s = s, w = w, epsilon = epsilon, W_p = W_p, maxit = maxit)
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
            new_theta[i,] <- linesearch_conic(theta = (new_theta[i-1,] + new_theta[i+1,])/2, s = s/2^k, w = w, epsilon = epsilon, W_p = W_p)
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
# plot(mesh)

region_conic <- function(theta_0, phi_l = 0, phi_r, eta, s = 1, epsilon, epsilon_tilde = epsilon, delta, tau = epsilon/1000, W_p, maxit = 1000, buff_boundary = FALSE, return_list = FALSE){
    w_up <- -log(delta/2)/eta
    w_lwr <- -log(1-(delta+tau)/2)/eta + W_p(theta_0 + tau * c(cos((phi_r + phi_l)/2), sin((phi_r+phi_l)/2)))
    theta <- pathfind_conic(phi_l = phi_l, phi_r = phi_r, s = s, epsilon = epsilon, w = w_up, epsilon_tilde =  epsilon_tilde, W_p = W_p, maxit = maxit, buff = tau)
    theta_small <- pathfind_conic(phi_l = phi_l, phi_r = phi_r, s = s, epsilon = epsilon/10, w = w_lwr, epsilon_tilde =  epsilon_tilde/100, W_p = W_p, maxit = maxit, buff = tau)
    theta_small <- theta_small[nrow(theta_small):1,]

    if(buff_boundary){
        theta <- cut_boundary(theta = theta, tau = tau)
    }

    boundary = rbind(theta, theta_small)


    if(return_list){
        return(list(inner = theta_small, outer = theta, boundary = boundary))
    } else{
        return(rbind(theta, theta_small))
    }
}


#' @noRd
# type = "h" or "v" -> horizontal or vertical
# direction = "positive" or "negative"
# W_p <- function(x){x[1]/(1-x[2])}
# linesearch_strip(c(0, 0.5), w = 0.6, epsilon = 1e-3, W_p = W_p)

linesearch_strip <- function(theta, type = "h", direction = "positive", s = 1, w, epsilon, W_p, maxit = 1000){
    if(type == "h"){
        v <- c(1, 0)
    } else{
        v <- c(0, 1)
    }
    direction <- ifelse(direction == "positive", 1, -1)
    it <- 0
    it2 <- 0
    while(abs(W_p(theta)-w) > epsilon^2 && it2 < maxit){
        while(W_p(theta) > w && it < maxit){
            theta <- theta - direction * s*v
            it <- it + 1
        }
        while((W_p(theta)<w) && it < maxit){
            theta <- theta + direction *s*v
            it <- it + 1
        }
        it <- 0
        it2 <- it2 + 1
        s <- s/2
    }
    if((it == maxit) || (it2 == maxit)){
        stop("Maximum number of iterations reached. Consider increasing 'maxit' to try again. Note that it is possible that W_p does not satisfy the required assumptions.")
    }
    return(theta)
}

#' @noRd 
# function to obtain the approximate level curve
# W_p <- function(x){x[1]/(1-x[2])}
# path <- pathfind_strip(start = 0, lower_bnd = 0, upper_bnd = 1, s = 1, epsilon = 1e-2, w = 0.5, epsilon_tilde = 1e-2, W_p = W_p)

pathfind_strip <- function(start, lower_bnd, upper_bnd, type = "h", direction = "positive", s, epsilon, w, epsilon_tilde, W_p, maxit = 1000, buff = 1e-6){
    direction_number <- ifelse(direction == "positive",1,-1)
    if(type == "h"){
        v1 <- c(start + direction_number * s, lower_bnd + buff)
        v2 <- c(start + direction_number * s, upper_bnd - buff)
    } else{
        v1 <- c(lower_bnd + buff, start + direction_number * s)
        v2 <- c(upper_bnd - buff, start + direction_number * s)
    }
    theta_0 <- linesearch_strip(theta = v1, type = type, direction = direction, s = s, w = w, epsilon = epsilon, W_p = W_p, maxit = maxit)
    theta_1 <- linesearch_strip(theta = v2, type = type, direction = direction, s = s, w = w, epsilon = epsilon, W_p = W_p, maxit = maxit)
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
            new_theta[i,] <- linesearch_strip(theta = (new_theta[i-1,] + new_theta[i+1,])/2, direction = direction, s = s/2^k, w = w, epsilon = epsilon, W_p = W_p)
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
# plot(mesh)
# start -> if type = "h", coordinate x of the base model, if type = "v", coordinate y of the base model
# lower_bnd/upper_bnd -> if type = "h", lower and upper bounds for the y coordinate, if type = "v", lower and upper bounds for the x coordinate

# # Examples:
# library(fmesher)
# W_p <- function(x){abs(x[1])/(1-x[2])}
# # example positive direction:
# reg <- region_strip(c(0,1), lower_bnd = 0, upper_bnd = 1, eta = 1, direction = "positive", s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p, tau = 0.01)
# mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
# plot(mesh)
# # example negative direction
# reg <- region_strip(c(0,1), lower_bnd = 0, upper_bnd = 1, eta = 1, direction = "negative", s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p, tau = 0.01)
# mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
# plot(mesh)
# # example both directions
# reg <- region_strip(c(0,1), lower_bnd = 0, upper_bnd = 1, eta = 1, direction = "both", s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p, tau = 0.01)
# mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
# plot(mesh)

region_strip <- function(theta_0, lower_bnd = 0, upper_bnd, eta, type = "h", direction = "positive", s = 1, epsilon, epsilon_tilde = epsilon, delta, tau = epsilon/1000, W_p, maxit = 1000, buff_boundary = FALSE, return_list = FALSE){
    if(type == "h"){
        v <- c(theta_0[1], (lower_bnd + upper_bnd)/2 + tau)
        start <- theta_0[1]
    } else{
        v <- c((lower_bnd + upper_bnd)/2 + tau, theta_0[2])
        start <- theta_0[2]
    }
    
    w_up <- -log(delta/2)/eta
    w_lwr <- -log(1-(delta+tau)/2)/eta + W_p(theta_0 + v)

    if(direction != "both"){    
        theta <- pathfind_strip(start = start, lower_bnd = lower_bnd, upper_bnd = upper_bnd, type = type, direction = direction, s = s, epsilon = epsilon, w = w_up, epsilon_tilde =  epsilon_tilde, W_p = W_p, maxit = maxit, buff = tau/10)
        theta_small <- pathfind_strip(start = start, lower_bnd = lower_bnd, upper_bnd = upper_bnd, type = type, direction = direction, s = s, epsilon = epsilon, w = w_lwr, epsilon_tilde =  epsilon_tilde, W_p = W_p, maxit = maxit, buff = tau/10)
    } else{
        theta_pos <- pathfind_strip(start = start, lower_bnd = lower_bnd, upper_bnd = upper_bnd, type = type, direction = "positive", s = s, epsilon = epsilon, w = w_up, epsilon_tilde =  epsilon_tilde, W_p = W_p, maxit = maxit, buff = tau/10)
        theta_small_pos <- pathfind_strip(start = start, lower_bnd = lower_bnd, upper_bnd = upper_bnd, type = type, direction = "positive", s = s, epsilon = epsilon/10, w = w_lwr, epsilon_tilde =  epsilon_tilde/10, W_p = W_p, maxit = maxit, buff = tau/10)
        theta_small_neg <- pathfind_strip(start = start, lower_bnd = lower_bnd, upper_bnd = upper_bnd, type = type, direction = "negative", s = s, epsilon = epsilon/10, w = w_lwr, epsilon_tilde =  epsilon_tilde/10, W_p = W_p, maxit = maxit, buff = tau/10)
        theta_neg <- pathfind_strip(start = start, lower_bnd = lower_bnd, upper_bnd = upper_bnd, type = type, direction = "negative", s = s, epsilon = epsilon, w = w_up, epsilon_tilde =  epsilon_tilde, W_p = W_p, maxit = maxit, buff = tau/10)
    }
    
    if(buff_boundary){
        if(direction != "both"){
            theta <- cut_boundary(theta = theta, tau = tau)
            theta_small <- cut_boundary(theta = theta_small, tau = tau)
        } else{
            theta_pos <- cut_boundary(theta = theta_pos, tau = tau)
            theta_small_pos <- cut_boundary(theta = theta_small_pos, tau = tau)
            theta_neg <- cut_boundary(theta = theta_neg, tau = tau)
            theta_small_neg <- cut_boundary(theta = theta_small_neg, tau = tau)
        }
    }

    if(direction == "negative"){
        theta <- theta[nrow(theta):1,]
    } else if(direction == "positive"){
        theta_small <- theta_small[nrow(theta_small):1,]
    } else{
        theta_small_pos <- theta_small_pos[nrow(theta_small_pos):1,]
        theta_neg <- theta_neg[nrow(theta_neg):1,]
    }

    if(direction == "both"){
        boundary = rbind(theta_pos, theta_small_pos, theta_small_neg, theta_neg)
    } else{
        boundary = rbind(theta, theta_small)
    }



    if(return_list){
        if(direction!="both"){
            return(list(inner = theta_small, outer = theta, boundary = boundary))
        } else{
            return(list(inner_pos = theta_small_pos, outer_pos = theta_pos, inner_neg = theta_small_neg, outer_neg = theta_neg, boundary = boundary))
        }
    } else{
        return(boundary)
    }
}


