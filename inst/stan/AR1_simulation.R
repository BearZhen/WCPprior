

rm(list = ls())
library(rstan)
WCP_coverage_count_phi = 0

PC_coverage_count_phi = 0

WCP_coverage_length_phi = c()

WCP_function_coverage_length_phi = c()

for (k in 1:100){
  print(k)

# simulate AR1 process
# n is the length of the process
n = 100
sigma = 0.1
phi = 0.8
sim_data = numeric(n)
sim_data[1] = rnorm(1,mean = 0, sd = sigma)
for (i in 2:n){
  sim_data[i] = phi*sim_data[i-1] + rnorm(1, mean = 0, sd = sigma*sqrt(1-phi^2) )
}
data <- list(N = 100L, 
             y = sim_data
)

fit_WCP_function <- stan(file = 'AR1_WCP_function_version.stan', data = data, iter = 1000)
WCP_function_post_phi <- extract(fit_WCP_function, pars = "phi")
WCP_function_hdi_phi = bayestestR::hdi(WCP_function_post_phi$phi, ci = 0.9)
if(phi >= WCP_function_hdi_phi$CI_low && phi <= WCP_function_hdi_phi$CI_high ){
  WCP_function_coverage_count_phi = WCP_function_coverage_count_phi + 1
  WCP_function_coverage_length_phi = c(WCP_function_coverage_length_phi, WCP_function_hdi_phi$CI_high - WCP_function_hdi_phi$CI_low)
}

fit_WCP <- stan(file = 'AR1_WCP.stan', data = data, iter = 1000)
WCP_post_phi <- extract(fit_WCP, pars = "phi")
WCP_hdi_phi = bayestestR::hdi(WCP_post_phi$phi, ci = 0.9)
if(phi>=WCP_hdi_phi$CI_low && phi<=WCP_hdi_phi$CI_high ){
  WCP_coverage_count_phi = WCP_coverage_count_phi + 1
  WCP_coverage_length_phi = c(WCP_coverage_length_phi, WCP_hdi_phi$CI_high - WCP_hdi_phi$CI_low)
}

}

WCP_average_cl_phi = mean(WCP_coverage_length_phi)
WCP_function_average_cl_phi = mean(WCP_function_coverage_length_phi)

#library(shinystan)
#launch_shinystan(fit_WCP)
