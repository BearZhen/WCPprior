functions{
#include /include/WCP_analytic.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower = -1, upper = 1> phi;
  //real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (n in 2:N)
    y[n] ~ normal(phi * y[n-1], 0.1);
  //eta = -3.932725e-11;
  phi ~ WCP_1D_AR1_analytic(N, 2.172101, 0.1);
}