
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
  target += log( abs(-0.1/sqrt(2)*((N*phi^N-1+phi^N-N*phi)*(1-phi)/(sqrt(N*(1-phi^2)-2*phi*(1-phi^N))) + sqrt(N*(1-phi^2)-2*phi*(1-phi^N)) )/((1-phi)^2*sqrt(N-sqrt(N*(1-phi^2)-2*phi*(1-phi^N))/(1-phi)))   ) );
  target += log( 2.172101*exp(-2.172101*sqrt(2)*0.1*sqrt(N - sqrt(N*(1-phi^2)-2*phi*(1-phi^N))/(1-phi)))/(1 - exp(-2.172101*sqrt(2)*0.1*sqrt(N-sqrt((1-(-1)^N)/2)) )) );
}

