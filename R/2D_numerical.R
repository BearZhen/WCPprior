
#' The 2d numerical approximation of WCP prior
#'
#' @param cutoff Cutoff parameter indicating probability tail mass for the Wasserstein distance.
#' @param W_func A function that returns the Wasserstein distance given theta1 and theta2. 
#' @param mesh_width Mesh width.
#' @param alpha Power of mesh width for a finer mesh.
#' @param tau tau parameter, indicating how much bounded away from boundary.
#' @param region A list that specify region information. For example, for strip regions (regions of the type IxJ, where I = (a,b) and J = (c,d) are intervals, with one of them being unbounded), one can set `region = list(type = "strip", base_theta = c(a,c), corners = c(a,b,c,d))`, where some of the values `a`,`b`,`c` or `d` can be `-Inf` or `Inf`. If the type is conic, the list must be set, for example, as list(type = 'conic', lower_angle = 0, upper_angle = pi, base_theta = c(0,0)), where `lower_angle` and `upper_angle` contains the smaller and larger angles that determine the cone.
#' @param eta User specified parameter of the WCP prior.
#' @param parallel A logic value that indicating whether the user wants to run the function with multiple cpu.
#' @param lc_multiplier Multiplier determines number of level curves
#' @param NumCores Number of cores to run the function.
#' @param visual_mesh A logic value to determin whether to plot the mesh.
#'
#' @return A list of density locations, densities evaluated on those locations and other utilities that can help access accuracy.
#' @export
#'
WCP_2D_Numerical_Density = function(W_func,
                          eta,
                          mesh_width, 
                          alpha = 1,
                          tau = mesh_width/1000,
                          cutoff = 0.01, 
                          region = list(type = 'conic', lower_angle = 0, upper_angle = pi, base_theta = c(0,0)),
                          #region = list(type = 'strip', corner = c(l1,u1,l2,u2), base_theta = c(0,0))
                          lc_multiplier = 20,
                          parallel = FALSE,
                          NumCores = parallel::detectCores()-1,
                          visual_mesh = FALSE) {
  
  ######################################## mesh generation ######################################
  ## mesh for interpolating the Wasserstein distance, partial derivatives of the Wasserstein distance
  if (region$type == 'conic'){
    # boundary path

    boundary_path = region_conic(theta_0 = region$base_theta, phi_l = region$lower_angle, phi_r = region$upper_angle, eta = eta, s = 1, epsilon = mesh_width, delta = cutoff, W_p = W_func)
    level_curve_type = 'LL'
    # lower bound of the second coordinate
    
  } else if (region$type == 'strip'){
    # determine horozontal or verticle and direction
    if (  !is.finite(region$corner[1]) | !is.finite(region$corner[2])   ){
      h_or_v = 'h'
      L = region$corner[3]
      U = region$corner[4]
      if ( is.finite(region$corner[1]) ){
        direction = 'positive'
        level_curve_type = 'L2U2'
      } else if ( !is.finite(region$corner[1]) & is.finite(region$corner[2]) ) {
        direction = 'negative'
        level_curve_type = 'L2U2'
      } else{
        direction = 'both'
        level_curve_type = 'L2L2'
      }
    } else {
      h_or_v = 'v'
      L = region$corner[1]
      U = region$corner[2]
      if ( is.finite(region$corner[3]) ){
        direction = 'positive'
        level_curve_type = 'L1U1'
      } else if ( !is.finite(region$corner[1]) & is.finite(region$corner[2]) ) {
        direction = 'negative'
        level_curve_type = 'L1U1'
      } else{
        
        direction = 'both'
        level_curve_type = 'L1L1'
      }
      
    }
    
    boundary_path = region_strip(theta_0 = region$base_theta, lower_bnd = region$corner[3], upper_bnd = region$corner[4], eta = eta, type = h_or_v, direction = direction, s = 1, epsilon = mesh_width, delta = cutoff, W_p = W_func, tau = tau)
    
 
    
  } else {
    stop("Wronly specifying the region or the type of region is not yet developed. Please contact our development team if you need.")
  }
  
  mesh = fmesher::fm_mesh_2d(boundary = fmesher::fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width)
  if (visual_mesh){
    plot(mesh)
  }
  
  # weights on mesh
  weights = numeric()
  for (i in 1:dim(mesh$loc)[1]){
    weights[i] = W_func(c(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2])))
  }
  
  
  # the default weights_fine and mesh_finer since the default alpha is 1.
  mesh_finer = mesh
  weights_fine = weights
  
  ## mesh for computing partial arc length values if alpha is not 1. Weights of this mesh will also be updated.
  if (alpha != 1){
    mesh_finer = fmesher::fm_mesh_2d(boundary = fmesher::fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width^alpha)
    weights_fine = numeric()
    for (i in 1:dim(mesh_finer$loc)[1]){
      weights_fine[i] = W_func( c(as.numeric(mesh_finer$loc[i,1]), as.numeric(mesh_finer$loc[i,2])) ) 
    }
  }
  
  
  ############# obtain level curves and filter the non-complete ones ###########################
  # get lower and upper bound of the Wasserstein distance
  W_lower_bound = min(weights_fine)
  W_upper_bound = max(weights_fine)

  # a vector containing partial arc length value on each discrete points representing level curves
  parc = numeric()
  # a vector containing total arc length value on each discrete point representing level curves
  tarc = numeric()
  # coordinates of discrete points representing level curves, these are the locations where we will evaluate prior density.
  density_location = numeric()
  
  if (parallel == TRUE){
  doParallel::registerDoParallel(cl <- parallel::makeCluster(NumCores))
  parallel::clusterEvalQ(cl, library("excursions"))
  parallel::clusterExport(cl, "weights_fine",
                          envir = environment())
  parallel::clusterExport(cl, "mesh",
                          envir = environment())
  parallel::clusterExport(cl, "compute_partial_arc_lengths",
                          envir = environment())
  
  results <- foreach::foreach (W = seq(from = W_lower_bound, to = W_upper_bound, length = lc_multiplier*ceiling(W_upper_bound/(mesh_width^alpha)))) %dopar% {
    levelcurve = excursions::tricontourmap(mesh_finer, z = weights_fine,
                               levels = W)$contour 
    # skip if there is no such level curve
    temp = try(sp::coordinates(levelcurve)[[1]][[1]], silent = FALSE)
    if ('try-error' %in% class(temp)){
      return(NULL)
    } else{ # otherwise, check the type pf level curve and update partial and total arc length
      
      # obtain the coordinates of all the discrete points of the level curve
      line_coord = sp::coordinates(levelcurve)[[1]][[1]]
      
      # if the second coordinate of level curves starts from L2 and ends at L2, like the Gaussian case
      if ( level_curve_type == 'LL'  ){
        # drop the non-complete curve
        # line_coord[1,2] - L2 > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - L2 > 2*mesh_width^alpha
        if (region$lower_angle < pi/2){
          L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] )
        } else{
          L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] ) + pi
        }
        if (region$upper_angle < pi/2){
          U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] )
        } else{
          U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] ) + pi
        }
        
        if ( abs(L_angle_line_coord - region$lower_angle) > 0.1*(region$upper_angle - region$lower_angle) | abs(U_angle_line_coord - region$upper_angle) > 0.1*(region$upper_angle - region$lower_angle) ) {
          return(NULL)
          
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) 
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
        
      } else if (level_curve_type == 'L2U2' & direction == 'positive' ){ # if the second coordinate of level curves starts from L2 and ends at U2, like the generalized Pareto case
        
        if (line_coord[1,2] - L > 5*mesh_width^alpha | U - line_coord[length(line_coord[,1]),2]  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + U - max(boundary_path[,2])
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
        
      } else if (level_curve_type == 'L2U2' & direction == 'negative'  ){
        if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(U - line_coord[1,2])  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + U - max(boundary_path[,2])
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) +  min(boundary_path[,2])
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
      } else if (level_curve_type == 'L2L2'){
        if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(L - line_coord[1,2])  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2])
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
      } else if (level_curve_type == 'L1U1' & direction == 'positive'){
        if ( abs(line_coord[1,1] - U) > 5*mesh_width^alpha | abs(L - line_coord[length(line_coord[,1]),1])  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + abs( U - max(boundary_path[,1]) )
          # update total arc length, compensated by the lift
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + abs( min(boundary_path[,1]) - L)
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }
      } else{
        stop("This case is under construction, please contact our development team.")
      }
      
    }
    
  }
  # obtain density locations, partial arc length and total arc length on density locations
  complete_results <- do.call(rbind, results)
  density_location <- complete_results[,1:2]
  parc <- complete_results[,3]
  tarc <- complete_results[,4]
  } else {
    
    for (W in seq(from = W_lower_bound, to = W_upper_bound, length = lc_multiplier*ceiling(W_upper_bound/(mesh_width^alpha)))){
    levelcurve = excursions::tricontourmap(mesh_finer, z = weights_fine,
                               levels = W)$contour 
    # skip if there is no such level curve
    temp = try(sp::coordinates(levelcurve)[[1]][[1]], silent = FALSE)
    if ('try-error' %in% class(temp)){
      return(NULL)
    } else{ # otherwise, check the type pf level curve and update partial and total arc length
      # obtain the coordinates of all the discrete points of the level curve
      line_coord = sp::coordinates(levelcurve)[[1]][[1]]
      
      # if the second coordinate of level curves starts from L2 and ends at L2, like the Gaussian case
      if ( level_curve_type == 'LL'  ){
        # drop the non-complete curve
        # line_coord[1,2] - L2 > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - L2 > 2*mesh_width^alpha
        if (region$lower_angle < pi/2){
          L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] )
        } else{
          L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] ) + pi
        }
        if (region$upper_angle < pi/2){
          U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] )
        } else{
          U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] ) + pi
        }
      
        if ( abs(L_angle_line_coord - region$lower_angle) > 0.1*(region$upper_angle - region$lower_angle) | abs(U_angle_line_coord - region$upper_angle) > 0.1*(region$upper_angle - region$lower_angle) ) {
          return(NULL)
          
        } else{
          # update density location
          density_location = rbind(density_location, line_coord)
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          parc = c(parc, levelcurve_parc)
          # update total arc length, compensated by the lift
          tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]))
        }
         
        
        
      } else if (level_curve_type == 'L2U2' & direction == 'positive' ){ # if the second coordinate of level curves starts from L2 and ends at U2, like the generalized Pareto case
        
        if (line_coord[1,2] - L > 5*mesh_width^alpha | U - line_coord[length(line_coord[,1]),2]  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update density location
          density_location = rbind(density_location, line_coord)
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          parc = c(parc, levelcurve_parc)
          # update total arc length, compensated by the lift
          tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc))  + U - max(boundary_path[,2]) )
        }
        
      } else if (level_curve_type == 'L2U2' & direction == 'negative'  ){
        if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(U - line_coord[1,2])  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update density location
          density_location = rbind(density_location, line_coord)
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + U - max(boundary_path[,2])
          parc = c(parc, levelcurve_parc)
          # update total arc length, compensated by the lift
          tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) +  min(boundary_path[,2]) )

        }
      } else if (level_curve_type == 'L2L2'){
        if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(L - line_coord[1,2])  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update density location
          density_location = rbind(density_location, line_coord)
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          parc = c(parc, levelcurve_parc)
          # update total arc length, compensated by the lift
          tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) )
 
        }
      } else if (level_curve_type == 'L1U1' & direction == 'positive'){
        if ( abs(line_coord[1,1] - U) > 5*mesh_width^alpha | abs(L - line_coord[length(line_coord[,1]),1])  > 5*mesh_width^alpha){
          return(NULL)
        } else{
          # update density location
          density_location = rbind(density_location, line_coord)
          # update partial arc length, compensated by the lift
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + abs( U - max(boundary_path[,1]) )
          parc = c(parc, levelcurve_parc)
          # update total arc length, compensated by the lift
          tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + abs( min(boundary_path[,1]) - L) )
        }
      } else{
        stop("This case is under construction, please contact our development team.")
      }
      
     # 
    }
    
    
    #
    }

  }






  ## check if every triangular element contains at least one point from level curves
  PD_util_P = fmesher::fm_basis(mesh, density_location, derivatives = TRUE)
  A = PD_util_P$A
  # find the mesh nodes whose element does not contain any level curve points
  idx_problem <- which(Matrix::colSums(A)==0)
  jitter = 0
  while(length(idx_problem)>0){
    new_loc <- mesh$loc[idx_problem, , drop=FALSE]
    new_W <- W_func( c(new_loc[,1], new_loc[,2])) + jitter
    if (parallel == TRUE){
    results_additional <- foreach::foreach (W = new_W) %dopar% {
      levelcurve = excursions::tricontourmap(mesh_finer, z = weights_fine,
                                 levels = W)$contour 
      # skip if there is no such level curve
      temp = try(sp::coordinates(levelcurve)[[1]][[1]], silent = FALSE)
      if ('try-error' %in% class(temp)){
        return(NULL)
      } else{ 
        # otherwise, check the type pf level curve and update partial and total arc length
        
        # obtain the coordinates of all the discrete points of the level curve
        line_coord = sp::coordinates(levelcurve)[[1]][[1]]
        
        # if the second coordinate of level curves starts from L2 and ends at L2, like the Gaussian case
        if ( level_curve_type == 'LL'  ){
          # drop the non-complete curve
          # line_coord[1,2] - L2 > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - L2 > 2*mesh_width^alpha
          if (region$lower_angle < pi/2){
            L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] )
          } else{
            L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] ) + pi
          }
          if (region$upper_angle < pi/2){
            U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] )
          } else{
            U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] ) + pi
          }
          
          if ( abs(L_angle_line_coord - region$lower_angle) > 0.1*(region$upper_angle - region$lower_angle) | abs(U_angle_line_coord - region$upper_angle) > 0.1*(region$upper_angle - region$lower_angle) ) {
            return(NULL)
            
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) 
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
          
        } else if (level_curve_type == 'L2U2' & direction == 'positive' ){ # if the second coordinate of level curves starts from L2 and ends at U2, like the generalized Pareto case
          
          if (line_coord[1,2] - L > 5*mesh_width^alpha | U - line_coord[length(line_coord[,1]),2]  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + U - max(boundary_path[,2])
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
          
        } else if (level_curve_type == 'L2U2' & direction == 'negative'  ){
          if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(U - line_coord[1,2])  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + U - max(boundary_path[,2])
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) +  min(boundary_path[,2])
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
        } else if (level_curve_type == 'L2L2'){
          if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(L - line_coord[1,2])  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2])
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
        } else if (level_curve_type == 'L1U1' & direction == 'positive'){
          if ( abs(line_coord[1,1] - U) > 5*mesh_width^alpha | abs(L - line_coord[length(line_coord[,1]),1])  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + abs( U - max(boundary_path[,1]) )
            # update total arc length, compensated by the lift
            tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + abs( min(boundary_path[,1]) - L)
            return(cbind(line_coord, levelcurve_parc,tarc,W))
          }
        } else{
          stop("This case is under construction, please contact our development team.")
        }
        
        
      }
      
      
    }
    complete_results <-  do.call(rbind, results_additional)
    density_location <- rbind(density_location, complete_results[,1:2])
    parc <- c(parc, complete_results[,3])
    tarc <- c(tarc,complete_results[,4])

    } else {
      for (W in new_W){
      levelcurve = excursions::tricontourmap(mesh_finer, z = weights_fine,
                                 levels = W)$contour 
      # skip if there is no such level curve
      temp = try(sp::coordinates(levelcurve)[[1]][[1]], silent = FALSE)
      if ('try-error' %in% class(temp)){
        return(NULL)
      } else{
        
        # otherwise, check the type pf level curve and update partial and total arc length
        # obtain the coordinates of all the discrete points of the level curve
        line_coord = sp::coordinates(levelcurve)[[1]][[1]]
        
        # if the second coordinate of level curves starts from L2 and ends at L2, like the Gaussian case
        if ( level_curve_type == 'LL'  ){
          # drop the non-complete curve
          # line_coord[1,2] - L2 > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - L2 > 2*mesh_width^alpha
          if (region$lower_angle < pi/2){
            L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] )
          } else{
            L_angle_line_coord = atan( line_coord[1,2]/line_coord[1,1] ) + pi
          }
          if (region$upper_angle < pi/2){
            U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] )
          } else{
            U_angle_line_coord = atan( line_coord[length(line_coord[,1]),2]/line_coord[length(line_coord[,1]),1] ) + pi
          }
          
          if ( abs(L_angle_line_coord - region$lower_angle) > 0.1*(region$upper_angle - region$lower_angle) | abs(U_angle_line_coord - region$upper_angle) > 0.1*(region$upper_angle - region$lower_angle) ) {
            return(NULL)
            
          } else{
            # update density location
            density_location = rbind(density_location, line_coord)
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            parc = c(parc, levelcurve_parc)
            # update total arc length, compensated by the lift
            tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]))
          }
          
          
          
        } else if (level_curve_type == 'L2U2' & direction == 'positive' ){ # if the second coordinate of level curves starts from L2 and ends at U2, like the generalized Pareto case
          
          if (line_coord[1,2] - L > 5*mesh_width^alpha | U - line_coord[length(line_coord[,1]),2]  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update density location
            density_location = rbind(density_location, line_coord)
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            parc = c(parc, levelcurve_parc)
            # update total arc length, compensated by the lift
            tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + U - max(boundary_path[,2]) )
                     
          }
          
        } else if (level_curve_type == 'L2U2' & direction == 'negative'  ){
          if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(U - line_coord[1,2])  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update density location
            density_location = rbind(density_location, line_coord)
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + U - max(boundary_path[,2])
            parc = c(parc, levelcurve_parc)
            # update total arc length, compensated by the lift
            tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) +  min(boundary_path[,2]) )
                     
          }
        } else if (level_curve_type == 'L2L2'){
          if ( abs(line_coord[length(line_coord[,1]),2] - L) > 5*mesh_width^alpha | abs(L - line_coord[1,2])  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update density location
            density_location = rbind(density_location, line_coord)
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
            parc = c(parc, levelcurve_parc)
            # update total arc length, compensated by the lift
            tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) )
                     
          }
        } else if (level_curve_type == 'L1U1' & direction == 'positive'){
          if ( abs(line_coord[1,1] - U) > 5*mesh_width^alpha | abs(L - line_coord[length(line_coord[,1]),1])  > 5*mesh_width^alpha){
            return(NULL)
          } else{
            # update density location
            density_location = rbind(density_location, line_coord)
            # update partial arc length, compensated by the lift
            levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + abs( U - max(boundary_path[,1]) )
            parc = c(parc, levelcurve_parc)
            # update total arc length, compensated by the lift
            tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + abs( min(boundary_path[,1]) - L) )
          }
        } else{
          stop("This case is under construction, please contact our development team.")
        }
        
        
        
        #  
      }
      }
    }
    
    PD_util_P = fmesher::fm_basis(mesh, density_location, derivatives = TRUE)
    A = PD_util_P$A
    idx_problem <- which(colSums(A)==0)
    jitter <- abs(stats::rnorm(1, mean = mesh_width/10, sd = mesh_width/10))
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
  PD_util_W = fmesher::fm_basis(mesh, density_location, derivatives = TRUE)
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
    W_distance[i] = W_func( c(as.numeric(density_location[i,1]), as.numeric(density_location[i,2])))
  }
  # the approximated density on level curve points
  approx_WCP_density = eta * approx_detJ_abs * exp(-eta * W_distance)/tarc
  mass_lump_C = fmesher::fm_fem(mesh, order = 1)$c0
  
  result = list(density_location, approx_WCP_density, mass_lump_C, A_new, W_distance, parc, tarc, approx_detJ_abs)
  return(result)
  
}

