functions {
#include ../inst/stan/include/WCP_analytic.stan
}
data {
  int<lower=0> N; //number of data
  real y[N]; //iid data
}

parameters {
  real<lower = 0> sigma;
  real<lower = 0, upper = 1> xi;
}
model {
  target += WCP1_2D_GP_analytic(sigma, xi, 10);
  target += generalized_Pareto(y, sigma, xi);
}
