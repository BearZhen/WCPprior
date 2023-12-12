// analytic 1D WCP prior for phi of AR1 process
real WCP_1D_AR1_analytic_log(real x, int n, real eta, real sigma){
    real density;
    real ldensity;
    real f;
    real c;
    real J;
<<<<<<< Updated upstream
    // for (i in 1:num_elements(x)){
       f = sqrt(n*(1 - x^2) - 2*x*(1 - x^n));
       c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) );
       J = sigma/sqrt(2) * ( (n*x^n - 1 + x^n - n*x)*(1 - x) + f^2 )/( f*sqrt(n - f/(1-x))*(1 - x)^2 );
       density = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - x) )  )/(1 - exp(-eta*c));
    // }
    ldensity = log(density);
||||||| Stash base
    for (i in 1:num_elements(x)){
       f = sqrt(n*(1 - x[i]^2) - 2*x[i]*(1 - x[i]^n));
       c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) );
       J = sigma/sqrt(2) * ( (n*x[i]^n - 1 + x[i]^n - n*x[i])*(1 - x[i]) + f^2 )/( f*sqrt(n - f/(1-x[i]))*(1 - x[i])^2 );
       density[i] = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - x[i]) )  )/(1 - exp(-eta*c));
    }
    ldensity = sum(log(density));
=======

    f = sqrt(n*(1 - x^2) - 2*x*(1 - x^n));
    c = sigma * sqrt( 2*n - sqrt(2)*sqrt(1 - (-1)^n) );
    J = sigma/sqrt(2) * ( (n*x^n - 1 + x^n - n*x)*(1 - x) + f^2 )/( f*sqrt(n - f/(1-x))*(1 - x)^2 );
    density = abs(J)*eta*exp(-2*eta*sigma^2* (n - f/(1 - x) )  )/(1 - exp(-eta*c));
    ldensity = log(density);
>>>>>>> Stashed changes
    return ldensity;
}




