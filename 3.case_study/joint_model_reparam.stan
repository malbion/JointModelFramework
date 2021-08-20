/* NDDM & RIM joint model

code to run on all focal species at once

Bimler 2021
*/ 
  
  
data {
  int<lower=1> S;          // number of species (rows in model matrix)
  int<lower=1> N;          // number of observations
  int<lower=0> K;          // number of neighbours (columns in model matrix)
  int<lower=0> I;          // number of interactions observed
  
  int<lower=0> species_ID[N];   // species index for observations
  int<lower=0> perform[N];      // response variable 
  
  int<lower=0> istart[S];       // indices matching species to interactions
  int<lower=0> iend[S];
  int<lower=0> icol[I];
  int<lower=0> irow[I];
  
  matrix[N,K] X;         // neighbour abundances (the model matrix)

} 

parameters {
  
  vector[S] beta_i0;    // species-specific intercept 
  vector<lower=0>[S] disp_dev; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (seeds)
  // disp_dev = 1/sqrt(phi)

  vector[I] beta_ij;     // vector of interactions which have been observed
  
  // vector<lower=0>[1] response1; // first species-specific response parameter
  // // constrained to be positive for identifiability
  // vector[S-1] responseSm1; // other species-specific response parameters
  
  // abov didn't work so let's try forcing all responses to be positive
  vector<lower=0>[S] response; // species-specific competitive response parameter
  // >= 0 to avoid bimodality in response and effect 

  unit_vector[K] effect; // species-specific effect parameter

  real<lower=0> sigma;	// scale parameter for the logistic distribution (used to estimate re's)
} 

transformed parameters {
  
  // transformed parameters constructed from paramaters above
  // vector[S] response;        // combined vector of species-specific responses
  matrix[S, K] inter_mat;    // the community interaction matrix
  vector[I] re;              // interactions as calculated by the re model
  vector[N] mu;              // the linear predictor for seed production
  
  inter_mat = rep_matrix(0, S, K); // fill the community interaction matrix with 0 (instead of NA)
  
  // match observed interactions parameters to the correct position in the community matrix
  for(s in 1:S) {
    for(i in istart[s]:iend[s]) {
      inter_mat[irow[i], icol[i]] = beta_ij[i];
    }
  }
  
  // individual fitness model 
  for(n in 1:N) {
       mu[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], inter_mat[species_ID[n], ]));  
  }

  // // stitch together the response values
  // response = append_row(response1, responseSm1);
  
  // build a vector of interaction parameters based on the response effect model
  for (i in 1:I) {
    re[i] = response[irow[i]]*effect[icol[i]];
  }
} 

model {

  // priors
  beta_i0 ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
  response ~ normal(0, 1);   // defining a lower limit aboves truncates the normal distribution
  // response1 ~ normal(0, 1);   // constrained by parameter defintion to be positive
  // responseSm1 ~ normal(0,1);
  // effect ~ normal(0, 1);   // we can use the default prior for the unit_vector effect
  sigma ~ cauchy(0, 1);

  // seed production, i.e. F_i
  for(n in 1:N) {
    perform[n] ~ neg_binomial_2(mu[n], (disp_dev[species_ID[n]]^2)^(-1));
    // in our case study, seed production shows a better fit to a negative binomial than poisson   
    // distribution
  }

  // response-effect interactions
  for (i in 1:I) {
    target += logistic_lpdf(re[i] | beta_ij[i], sigma);
  }
  
}
