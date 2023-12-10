
#' The 2d numerical approximation of WCP prior
#'
#' @param boundary_path A sequence of points that represent boundary of computational domain.
#' @param W_func A function that returns the Wasserstein distance given theta1 and theta2. 
#' @param mesh_width Mesh width.
#' @param alpha Power of mesh width for a finer mesh.
#' @param W_lower_bound The lower bound of the Wasserstein distance.
#' @param W_upper_bound The upper bound of the Wasserstein distance.
#' @param eta User specified parameter of the WCP prior.
#' @param L2 Lower bound of theta 2.
#' @param U2 Upper bound of theta 2.
#' @param level_curve_type Type of level curve, the value should be either LL or LU
#' @param lc_multiplier Multiplier determines number of level curves
#' @param NumCores Number of cores to run the function.
#'
#' @return A list of density locations, densities evaluated on those locations and other utilities that can help access accuracy.
#' @export
#'


WCP_2D_Numerical_Density = function(boundary_path,
                          W_func,
                          mesh_width,
                          alpha = 1,
                          W_lower_bound,
                          W_upper_bound,
                          eta,
                          L2, U2,
                          level_curve_type,
                          lc_multiplier = 20,
                          NumCores){
  
  ######################################## mesh generation ######################################
  ## mesh for interpolating the Wasserstein distance, partial derivatives of the Wasserstein distance
  mesh = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width)
  plot(mesh)
  # weights on mesh
  weights = numeric()
  for (i in 1:dim(mesh$loc)[1]){
    weights[i] = W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2]))
  }
  # the default weights_fine and mesh_finer since the default alpha is 1.
  mesh_finer = mesh
  weights_fine = weights
  
  ## mesh for computing partial arc length values if alpha is not 1. Weights of this mesh will also be updated.
  if (alpha != 1){
    mesh_finer = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width^alpha)
    weights_fine = numeric()
    for (i in 1:dim(mesh_finer$loc)[1]){
      weights_fine[i] = W_func(as.numeric(mesh_finer$loc[i,1]), as.numeric(mesh_finer$loc[i,2]))
    }
  }
  
  ############# obtain level curves and filter the non-complete ones ###########################
  # a vector containing partial arc length value on each discrete points representing level curves
  parc = numeric()
  # a vector containing total arc length value on each discrete point representing level curves
  tarc = numeric()
  # coordinates of discrete points representing level curves, these are the locations where we will evaluate prior density.
  density_location = numeric()
  
  registerDoParallel(cl <- makeCluster(NumCores))
  parallel::clusterEvalQ(cl, library("excursions"))
  parallel::clusterExport(cl, "weights_fine",
                          envir = environment())
  parallel::clusterExport(cl, "mesh",
                          envir = environment())
  parallel::clusterExport(cl, "compute_partial_arc_lengths",
                          envir = environment())
  
  results <- foreach (W = seq(from = W_lower_bound, to = W_upper_bound, length = lc_multiplier*ceiling(W_upper_bound/(mesh_width^alpha)))) %dopar% {
    levelcurve = tricontourmap(mesh_finer, z = weights_fine,
                               levels = W)$contour 
    # skip if there is no such level curve
    temp = try(coordinates(levelcurve)[[1]][[1]], silent = FALSE)
    if ('try-error' %in% class(temp)){
      return(NULL)
    } else{ # otherwise, check the type pf level curve and update partial and total arc length
      
      # obtain the coordinates of all the discrete points of the level curve
      line_coord = coordinates(levelcurve)[[1]][[1]]
      
      # if the second coordinate of level curves starts from L2 and ends at L2, like the Gaussian case
      if (level_curve_type == 'LL'){
        # drop the non-complete curve
        if (line_coord[1,2] - L2 > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - L2 > 2*mesh_width^alpha){
          return(NULL)
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) 
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
        
      } else { # if the second coordinate of level curves starts from L2 and ends at U2, like the generalized Pareto case
        
        if (line_coord[1,2] - L2 > 5*mesh_width^alpha | U2 - line_coord[length(line_coord[,1]),2]  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + U2 - max(boundary_path[,2])
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
        
      }
      
    }
    
  }
  # obtain density locations, partial arc length and total arc length on density locations
  complete_results <- do.call(rbind, results)
  density_location <- complete_results[,1:2]
  parc <- complete_results[,3]
  tarc <- complete_results[,4]
  
  ## check if every triangular element contains at least one point from level curves
  PD_util_P = fm_basis(mesh, density_location, derivatives = TRUE)
  A = PD_util_P$A
  # find the mesh nodes whose element does not contain any level curve points
  idx_problem <- which(colSums(A)==0)
  print("length")
  print(length(idx_problem))
  jitter = 0
  while(length(idx_problem)>0){
    print("having holes")
    new_loc <- mesh$loc[idx_problem, , drop=FALSE]
    new_W <- W_func(new_loc[,1], new_loc[,2]) + jitter
    results_additional <- foreach (W = new_W) %dopar% {
      print(W)
      levelcurve = tricontourmap(mesh_finer, z = weights_fine,
                                 levels = W)$contour 
      # skip if there is no such level curve
      temp = try(coordinates(levelcurve)[[1]][[1]], silent = FALSE)
      if ('try-error' %in% class(temp)){
        return(NULL)
      } else{ # otherwise, check the type pf level curve and update partial and total arc length
        
        # obtain the coordinates of all the discrete points of the level curve
        line_coord = coordinates(levelcurve)[[1]][[1]]
        
        # if the second coordinate of level curves starts from L2 and ends at L2, like the Gaussian case
        if (level_curve_type == 'LL'){
          # drop the non-complete curve
          if (line_coord[1,2] - L2 > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - L2 > 2*mesh_width^alpha){
            return(NULL)
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) 
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
          
        } else { # if the second coordinate of level curves starts from L2 and ends at U2, like the generalized Pareto case
          if (line_coord[1,2] - L2 > 10*(mesh_width^alpha) | U2 - line_coord[length(line_coord[,1]),2]  > 10*(mesh_width^alpha) ){
            return(NULL)
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2])
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
        }
        
      }
      
      
    }
    complete_results <-  do.call(rbind, results_additional)
    density_location <- rbind(density_location, complete_results[,1:2])
    parc <- c(parc, complete_results[,3])
    tarc <- c(tarc,complete_results[,4])
    PD_util_P = fm_basis(mesh, density_location, derivatives = TRUE)
    A = PD_util_P$A
    idx_problem <- which(colSums(A)==0)
    jitter <- abs(rnorm(1, mean = mesh_width/10, sd = mesh_width/10))
    print("Adding more levels")
  }
  
  ########################### compute partial derivatives of partial arc length and the Wasserstein distance #########
  # interpolate partial arc lengths on mesh nodes from partial arc lengths on level curves points
  A_new <- solve(t(A)%*%A,t(A))
  interpolated_parc <- A_new%*%parc
  # compute partial derivatives of partial arc lengths on level curve points
  P_PD_x = PD_util_P$dx %*% matrix(interpolated_parc, ncol = 1)
  P_PD_x = as.vector(P_PD_x)
  P_PD_y = PD_util_P$dy %*% matrix(interpolated_parc, ncol = 1)
  P_PD_y = as.vector(P_PD_y)
  
  # compute partial derivatives of the Wasserstein distance on level curve points
  PD_util_W = fm_basis(mesh, density_location, derivatives = TRUE)
  W_PD_x = PD_util_W$dx %*% matrix(weights,ncol = 1)
  W_PD_x = as.vector(W_PD_x)
  W_PD_y = PD_util_W$dy %*% matrix(weights,ncol = 1)
  W_PD_y = as.vector(W_PD_y)
  
  ################## assemble the approximated density on level curve points ##########
  # compute absolute value of determinant of Jacobian on level curve points
  approx_detJ_abs = abs(W_PD_x * P_PD_y - W_PD_y * P_PD_x)
  # compute the Wasserstein distance on level curve points
  W_distance = numeric()
  for (i in 1:dim(density_location)[1]){
    W_distance[i] = W_func(as.numeric(density_location[i,1]), as.numeric(density_location[i,2]))
  }
  # the approximated density on level curve points
  approx_WCP_density = eta * approx_detJ_abs * exp(-eta * W_distance)/tarc
  mass_lump_C = fm_fem(mesh, order = 1)$c0
  
  result = list(density_location, approx_WCP_density, mass_lump_C, A_new, W_distance, parc, tarc, approx_detJ_abs)
  return(result)
  
}

