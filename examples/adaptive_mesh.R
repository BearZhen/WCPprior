Adaptive_2D_grid = function(  theta_1_init
                              , theta_2_init
                              , mesh
                              , FEM_weights
                              , ratio 
                              , tol.ratio = 0.7
                              #,W_lower
                              #,W_upper
                              ){
    
    # theta_1_init is a vector of values of theta1 of the initial grid
    # theta_2_init is a vector of values of theta2 of the initial grid
    
    grid_coord = expand.grid(theta_1_init, theta_2_init)
    grid_coord = as.matrix(grid_coord)
    # wasserstein distance on the grid
    A = inla.spde.make.A( mesh, loc = grid_coord ) 
    W_grid_init = A %*% weights
    
    for(i in 1:length(theta_2_init)){
        print('i')
        print(i)
        for (j in 1:length(theta_1_init)){
            #print('j')
            #print(j)
            if (j==length(theta_1_init)){
                break
            }
            index = (i-1)*length(theta_2_init) + j
            inserted_W = W_grid_init[index+1] 
            step_size = abs(grid_coord[index,1] - grid_coord[index+1,1]) 
            #print(abs(W_grid_init[index] - inserted_W )/W_grid_init[index])
            while ( as.numeric(abs(W_grid_init[index] - inserted_W )/W_grid_init[index]) > tol.ratio){
             step_size = step_size * ratio
             inserted_x_coord = grid_coord[index,1] + step_size 
             inserted_coord = t(as.matrix(c(inserted_x_coord, theta_2_init[i])))
             A_matrix = inla.spde.make.A( mesh, loc = inserted_coord ) 
             inserted_W = A_matrix %*% weights
             #print(abs(W_grid_init[index] - inserted_W )/W_grid_init[index])
            }
            theta_1_init = c(theta_1_init, inserted_x_coord)
            theta_1_init = sort(theta_1_init)
            print(length(theta_1_init))
        }
    }
return(theta_1_init)

}
theta_1_init = seq(from = 0.001, to = 0.3, by = 0.01)
theta_2_init = seq(from = 0, to = 0.999, by = 0.01)
adagrid = Adaptive_2D_grid(  theta_1_init = theta_1
                              , theta_2_init = theta_1
                              , mesh = mesh
                              , FEM_weights = weights
                              , ratio = 0.7
                              , tol.ratio = 0.2
                              )





a = c(1,2,3)
for (i in a){
    print(i)
    if(i == 1){
        a = c(a,1.5)
        a = sort(a)
    }
}
