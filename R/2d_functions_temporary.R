

WCP_2D_Density = function(boundary_path,
                          W_func,
                          mesh_width,
                          alpha){
  
  
  ######################################## mesh generation ######################################
  ## mesh for interpolating the Wasserstein distance, partial derivatives of the Wasserstein distance
  mesh = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width)
  # weights on mesh
  weights = numeric()
  for (i in 1:dim(mesh_coarser$loc)[1]){
    weights[i] = W_func(as.numeric(mesh_coarser$loc[i,1]), as.numeric(mesh_coarser$loc[i,2]))
  }
  
  ## mesh for computing partial arc length values
  if (alpha != 1){
    mesh_finer = fm_mesh_2d(boundary = fm_segm( boundary_path, is.bnd = TRUE), max.edge = mesh_width^alpha)
    weights_fine = numeric()
    for (i in 1:dim(mesh$loc)[1]){
      weights_fine[i] = W_func(as.numeric(mesh$loc[i,1]), as.numeric(mesh$loc[i,2]))
    }
  }
  

  
  ############# obtain level curves and filter the non-complete ones ###########################
  # a vector contianing partial arc length value on each discrete points representing level curves
  parc = numeric()
  # a vector contianing total arc length value on each discrete point representing level curves
  tarc = numeric()
  # coordinates of discrete points representing level curves, these are the locations where we will evaluate prior density.
  density_location = numeric()

  
  registerDoParallel(cl <- makeCluster(20))
  parallel::clusterEvalQ(cl, library("excursions"))
  parallel::clusterExport(cl, "weights_fine",
                          envir = environment())
  parallel::clusterExport(cl, "mesh",
                          envir = environment())
  parallel::clusterExport(cl, "compute_partial_arc_lengths",
                          envir = environment())
  
  results <- foreach (W = seq(from = W_lower_bound, to = W_upper_bound, length = 20*ceiling(W_upper_bound/(mesh_width^alpha)))) %dopar% {
    print(W)
    levelcurve = tricontourmap(mesh, z = weights_fine,
                               #tol = 1e-6,
                               levels = W)$contour #seq(from = W_lower_bound, to = W_upper_bound, length = ceiling(W_upper_bound/mesh_width)))$contour
    # skip if there is no such level curve
    temp = try(coordinates(levelcurve)[[1]][[1]], silent = FALSE)
    if ('try-error' %in% class(temp)){
      return(NULL)
    } else{
      # obtain the coordinates of all the discrete points of the level curve
      line_coord = coordinates(levelcurve)[[1]][[1]]
      # drop the non-complete curve, need to be adjust case by case
      if (line_coord[1,2] - theta2_low > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - theta2_low > 2*mesh_width^alpha){
        return(NULL)
      } else{
        # update density location
        # density_location = rbind(density_location,line_coord)
        # update 
        levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
        # (pi-acos(line_coord[,1]/(sqrt(line_coord[,1]^2 + line_coord[,2]^2))))*(sqrt(line_coord[,1]^2 + line_coord[,2]^2))
        #parc = c(parc, levelcurve_parc)
        # adjustment of partial arc length, need to be done case by case
        # tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) )
        tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) 
        return(cbind(line_coord, levelcurve_parc,tarc,W))
      }}
  }
  
  complete_results <- do.call(rbind, results)
  
  density_location <- complete_results[,1:2]
  parc <- complete_results[,3]
  tarc <- complete_results[,4]
  
  # build a coarser mesh
  # 
  # mean(abs((pi-acos(line_coord[,1]/(sqrt(line_coord[,1]^2 + line_coord[,2]^2))))*(sqrt(line_coord[,1]^2 + line_coord[,2]^2))-levelcurve_parc))
  
  print("Error parc before")
  print(mean(abs((pi-acos(density_location[,1]/(sqrt(density_location[,1]^2 + density_location[,2]^2))))*(sqrt(density_location[,1]^2 + density_location[,2]^2))-parc)))
  
  PD_util_P = fm_basis(mesh_coarser, density_location, derivatives = TRUE)
  A = PD_util_P$A
  # 176739  11180
  # 29036 11180
  idx_problem <- which(colSums(A)==0)
  jitter = 0
  while(length(idx_problem)>0){
    new_loc <- mesh_coarser$loc[idx_problem, , drop=FALSE]
    new_W <- W_func(new_loc[,1], new_loc[,2]) + jitter
    results_additional <- foreach (W = new_W) %dopar% {
      print(W)
      levelcurve = tricontourmap(mesh, z = weights_fine,
                                 #tol = 1e-6,
                                 levels = W)$contour #seq(from = W_lower_bound, to = W_upper_bound, length = ceiling(W_upper_bound/mesh_width)))$contour
      # skip if there is no such level curve
      temp = try(coordinates(levelcurve)[[1]][[1]], silent = FALSE)
      if ('try-error' %in% class(temp)){
        return(NULL)
      } else{
        # obtain the coordinates of all the discrete points of the level curve
        line_coord = coordinates(levelcurve)[[1]][[1]]
        # drop the non-complete curve, need to be adjust case by case
        if (line_coord[1,2] - theta2_low > 2*mesh_width^alpha | line_coord[length(line_coord[,1]),2] - theta2_low > 2*mesh_width^alpha){
          return(NULL)
        } else{
          # update density location
          # density_location = rbind(density_location,line_coord)
          # update 
          levelcurve_parc = compute_partial_arc_lengths(line_coord)[,3] + min(boundary_path[,2])
          # (pi-acos(line_coord[,1]/(sqrt(line_coord[,1]^2 + line_coord[,2]^2))))*(sqrt(line_coord[,1]^2 + line_coord[,2]^2))
          #parc = c(parc, levelcurve_parc)
          # adjustment of partial arc length, need to be done case by case
          # tarc = c(tarc, rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) )
          tarc = rep(max(levelcurve_parc), times = length(levelcurve_parc)) + min(boundary_path[,2]) 
          return(cbind(line_coord, levelcurve_parc,tarc,W))
        }}
    }
    complete_results <-  do.call(rbind, results_additional)
    density_location <- rbind(density_location, complete_results[,1:2])
    parc <- c(parc, complete_results[,3])
    tarc <- c(tarc,complete_results[,4])
    PD_util_P = fm_basis(mesh_coarser, density_location, derivatives = TRUE)
    A = PD_util_P$A
    idx_problem <- which(colSums(A)==0)
    jitter <- abs(rnorm(1, mean = mesh_width/10, sd = mesh_width/10))
    print("Adding more levels")
  }
  
  print("Finished!")
  
}