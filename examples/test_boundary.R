
library(fmesher)
W_p <- function(x){sqrt(x[1]^2+x[2]^2)}
reg <- region_conic(theta_0 = c(0,0), phi_l = 0, phi_r = pi, eta = 1, s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p)
mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
plot(mesh)



W_p <- function(x){abs(x[1])/(1-x[2])}
# example positive direction:
reg <- region_strip(c(0,0), lower_bnd = 0, upper_bnd = 1, eta = 20, direction = "positive", s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p, tau = 0.001)
mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
plot(mesh, axe = TRUE)
# example negative direction
reg <- region_strip(c(0,0), lower_bnd = 0, upper_bnd = 1, eta = 20, direction = "negative", s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p, tau = 0.01)
mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
plot(mesh)
# example both directions
reg <- region_strip(c(0,1), lower_bnd = 0, upper_bnd = 1, eta = 1, direction = "both", s = 1, epsilon = 0.1, delta = 0.01, W_p = W_p, tau = 0.01)
mesh = fm_mesh_2d(boundary = fm_segm( reg, is.bnd = TRUE), max.edge = 0.1, cutoff = 0.01)
plot(mesh)
