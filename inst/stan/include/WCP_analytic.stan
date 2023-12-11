// analytic 1D WCP prior for phi of AR1 process
functions{
  real WCP_1D_AR1_analytic_log(vector x, int n, real eta, real sigma){
    vector[num_elements(x)] density;
    real ldensity;
    real f;
    real c;
    real J;
    for (i in 1:num_elements(x)){
       f = sqrt(n*(1 - x[i]^2) - 2*x[i]*(1 - x[i]^n));
       c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) );
       J = sigma/sqrt(2) * ( (n*x[i]^n - 1 + x[i]^n - n*x[i])*(1 - x[i]) + f^2 )/( f*sqrt(n - f/(1-x[i]))*(1 - x[i])^2 );
       density[i] = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - x[i]) )  )/(1 - exp(-eta*c));
    }
    ldensity = sum(log(density));
    return ldensity;
  }
}



