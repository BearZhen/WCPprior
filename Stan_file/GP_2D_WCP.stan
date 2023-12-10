
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; //number of data
  real y[N]; //iid data
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0, upper = 1> xi;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += log( pow(sigma,2.0) + pow(sigma, 4.0)/pow(1-xi, 2.0) );
  target += -log( pow(sigma, 2.0)*(1-xi)*sqrt(pow(sigma,2.0)/pow(1-xi,2.0) + 1) );
  target += log( 10*exp(- 10*sigma/(1-xi))/sqrt(pow(sigma,2.0)/pow(1-xi,2.0) + 1) );
  for(i in 1:N){
    target += log( pow(1.0 + xi*y[i]/sigma, -1/xi-1)/sigma );
   }
}

