// analytic 1D WCP2 prior for phi of AR1 process
real WCP_1D_AR1_analytic_log(real x, int n, real eta, real sigma){
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
real WCP_1D_GPtail_analytic_log(real x, real eta){
    real density;
    real ldensity;
    
    density =  eta/(1 - x)^2 * exp( - eta * x/(1 - x) );
    ldensity = log(density);

    return ldensity;
}

// analytic 1D WCP2 for for precision parameter tau (1/variance) of Gaussian distribution
real WCP_1D_Gaussian_precision_analytic_log(real x, real eta){
    real density;
    real ldensity;
    
    density =  0.5 * x^(-1.5) * eta * exp( -eta * x^( -0.5 ) );
    ldensity = log(density);

    return ldensity;
}

// analytic 1D WCP2 for for mean parameter of Gaussian distribution
real WCP_1D_Gaussian_mean_analytic_log(real x, real eta){
    real density;
    real ldensity;
    
    density =  0.5 * eta * exp( - eta * abs(x) );
    ldensity = log(density);

    return ldensity;
}





