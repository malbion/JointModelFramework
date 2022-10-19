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
  
  real<lower=0> weight;    // weighting value controlling the average strength of interactions (must be positive)
  unit_vector[S] response; // species-specific effect parameter (r in the manuscript but without weight)

  unit_vector[T] effect; // species-specific effect parameter

} 

transformed parameters {
  
  // transformed parameters constructed from parameters above
  vector[N] mu;              // the linear predictor for perform (here seed production)
  
  matrix[S, T] ri_betaij;   // interaction matrix for the RI model
  
  // get RIM interactions
  ri_betaij = weight * response*effect';  // the apostrophe transposes the effect vector
  
  // response-impact model estimates all interactions
  for(n in 1:N) {
       mu[n] = exp(gamma_i[species_ID[n]] - dot_product(X[n], ri_betaij[species_ID[n], ]));  
  }
  
} 

model {

  // priors
  gamma_i ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
 
  weight ~ normal(0, 10); // constrained by parameter definition to be positive
  // no prior needed for response or effect as we can use the default prior for the unit_vector

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




