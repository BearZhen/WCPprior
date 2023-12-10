
functions{
  real WCP_1D_AR1_analytic_log(vector x, real n, real eta, real sigma){
    vector[num_elements(x)] density;
    real ldensity;
    real f;
    real c;
    real J;
    for (phi in 1:num_elements(x)){
       f = sqrt(n*(1 - phi^2) - 2*phi*(1 - phi^n))
       c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) )
       J = sigma/sqrt(2) * ( (n*phi^n - 1 + phi^n - n*phi)*(1 - phi) + f^2 )/( f*sqrt(n - f/(1-phi))*(1 - phi)^2 )
       density[i] = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - phi) )  )/(1 - exp(-eta*c))
    }
  }
}


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  // length of the process
  int<lower = 1> n;
  real<lower = 0> sigma;
  real eta;
}


model {
  y ~ WCP_1D_AR1_analytic( n, eta, sigma);
}

