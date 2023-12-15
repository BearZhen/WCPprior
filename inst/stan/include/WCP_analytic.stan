// analytic 1D WCP2 prior for phi of AR1 process
real WCP2_1D_AR1_analytic_log(real x, int n, real eta, real sigma){
    real density;
    real ldensity;
    real f;
    real c;
    real J;
    f = sqrt(n*(1 - x^2) - 2*x*(1 - x^n));
    c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) );
    J = sigma/sqrt(2) * ( (n*x^n - 1 + x^n - n*x)*(1 - x) + f^2 )/( f*sqrt(n - f/(1-x))*(1 - x)^2 );
    density = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - x) )  )/(1 - exp(-eta*c));
    ldensity = log(density);

    return ldensity;
}

// analytic 1D WCP1 for xi (tail index) of generalized Pareto distribution
real WCP1_1D_GPtail_analytic_log(real x, real eta){
    real density;
    real ldensity;
    
    density =  eta/(1 - x)^2 * exp( - eta * x/(1 - x) );
    ldensity = log(density);

    return ldensity;
}

// analytic 1D WCP2 for for precision parameter tau (1/variance) of Gaussian distribution
real WCP2_1D_Gaussian_precision_analytic_log(real x, real eta){
    real density;
    real ldensity;
    
    density =  0.5 * x^(-1.5) * eta * exp( -eta * x^( -0.5 ) );
    ldensity = log(density);

    return ldensity;
}

// analytic 1D WCP2 for for mean parameter of Gaussian distribution
real WCP2_1D_Gaussian_mean_analytic_log(real x, real eta){
    real density;
    real ldensity;
    
    density =  0.5 * eta * exp( - eta * abs(x) );
    ldensity = log(density);

    return ldensity;
}

//  2d analytic density WCP prior for mean and standard deviation of Gaussian distribution
real WCP2_2D_Gaussian_analytic_log(vector x, real eta, real base_m){
  // x[1] is mean parameter, x[2] is standard deviation parameter
  real density;
  real ldensity;
  density = eta/sqrt( (x[1] - base_m )^2 + x[2]^2) * exp(-eta * sqrt( (x[1] - base_m )^2 + x[2]^2))/pi();
  ldensity = log(density);
  return ldensity;
}


//  2d analytic density WCP prior for sigma and xi of generalized Pareto distribution
real WCP1_2D_GP_analytic_log(real sigma, real xi, real eta){
  // x[1] is sigma, x[2] is xi
  real density;
  real ldensity;
  density = eta/(1 - xi) * exp(-eta * sigma/(1 - xi) );
  ldensity = log(density);
  return ldensity;
}

// log-likelihood of generalized Pareto
real generalized_Pareto_log(vector x, real sigma, real xi){
  vector [num_elements(x)] prob;
  real lprob;
  for (i in 1:num_elements(x)){
    prob[i] = (sigma^(-1)) * (1 + xi*x[i]/sigma)^(-1/xi - 1);
  }
  lprob = sum(log(prob));
  return lprob;
}



