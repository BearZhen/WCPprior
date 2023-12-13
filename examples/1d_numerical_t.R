library(ggplot2)
library(fGarch)

WpdistanceGt = function(xi){
  # xi is 1/nu where nu is dof of t-distribution
  integrand = function(p,...){
    abs(qstd(p,...) - qnorm(p))^2
  }
  integral = integrate(Vectorize(integrand),0,1, nu = 1/xi,subdivisions = 10000L)$value
  return (integral^(1/2))
}

Xi = seq(from = 0.001,to = 0.499,by = 0.001)
Wpdistance = c()
xi = c()
for (i in 1:length(Xi)){
  temp = try(WpdistanceGt(Xi[i]), silent = FALSE)
  if ('try-error' %in% class(temp)){
    print(Xi[i])
    next
  }
  Wpdistance = c(Wpdistance, WpdistanceGt(Xi[i]) )
  xi = c(xi, Xi[i])
}

W2d = splinefun(xi, Wpdistance, method = "hyman")

alpha = 0.3
U = 0.1
eta_WCP = -log(alpha)/W2d(U)

result = WCP_1D_Numerical_Density(base_theta = 0,
                                   L = 0, 
                                   L_included = TRUE,
                                   U = 0.5,
                                   U_included = FALSE,
                                   W_func = W2d,
                                   eta = eta_WCP,
                                   cutoff1 = NULL,
                                   cutoff2 = NULL,
                                   mesh_width = 0.001)

plot(result[[1]], result[[2]],type="l")
