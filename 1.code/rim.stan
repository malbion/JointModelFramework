/* RIM

code to run on all focal species at once

Bimler 2022
*/ 
  
  
data {
  int<lower=1> S;          // number of species (elements) 
  int<lower=1> N;          // number of observations (rows in model matrix)
  int<lower=0> T;          // number of neighbours (columns in model matrix)
  
  int<lower=0> species_ID[N];   // index matching species to observations
  int<lower=0> perform[N];      // response variable 
  
  matrix[N,T] X;         // neighbour abundances (the model matrix)

} 

parameters {
  
  vector[S] gamma_i;    // species-specific intercept 
    
  vector<lower=0>[S] disp_dev; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
  vector<lower=0>[1] response1; // first species-specific response parameter
  // constrained to be positive for identifiability
  vector[S-1] responseSm1; // other species-specific response parameters

  unit_vector[T] effect; // species-specific effect parameter

} 

transformed parameters {
  
  // transformed parameters constructed from parameters above
  vector[S] response;        // combined vector of species-specific responses
  
  vector[N] mu;              // the linear predictor for perform (here seed production)
  
  matrix[S, T] ri_betaij;   // interaction matrix for the RI model
  
  
  // stitch together the response values
  response = append_row(response1, responseSm1);
  // get RIM interactions
  ri_betaij = response*effect';
  
  // response-impact model estimates all interactions
  for(n in 1:N) {
       mu[n] = exp(gamma_i[species_ID[n]] - dot_product(X[n], ri_betaij[species_ID[n], ]));  
  }
  
} 

model {

  // priors
  gamma_i ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
 
  response1 ~ normal(0, 1);   // constrained by parameter defintion to be positive
  responseSm1 ~ normal(0,1);
  // no prior needed for effect as we can use the default prior for the unit_vector

  
  for(n in 1:N) {
    perform[n] ~ neg_binomial_2(mu[n], (disp_dev[species_ID[n]]^2)^(-1));
    // in our case study, seed production shows a better fit to a negative binomial 
    // than poisson distribution
  }
  
}

generated quantities {

  vector[N] log_lik;
  
  for(n in 1:N) {
    log_lik[n] = neg_binomial_2_lpmf(perform[n] | mu[n], (disp_dev[species_ID[n]]^2)^(-1));
  }

}




