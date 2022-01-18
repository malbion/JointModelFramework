/* NDDM & RIM joint model

code to run on all focal species at once

Bimler 2021
*/ 
  
  
data {
  int<lower=1> S;          // number of species (elements) 
  int<lower=1> N;          // number of observations (rows in model matrix)
  int<lower=0> K;          // number of neighbours (columns in model matrix)
  int<lower=0> I;          // number of inferred interactions 
  
  int<lower=0> species_ID[N];   // index matching species to observations
  int<lower=0> perform[N];      // response variable 
  //int<lower=0> perform2[N];      // response variable 
    
  int<lower=0> istart[S];       // indices matching species to inferred interactions
  int<lower=0> iend[S];
  int<lower=0> icol[I];
  int<lower=0> irow[I];
  
  matrix[N,K] X;         // neighbour abundances (the model matrix)
  // matrix[S,K] Q;	  // matrix of inferrable interactions

} 

parameters {
  
  vector[S] beta_i0;    // species-specific intercept 
    
  vector<lower=0>[S] disp_dev; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
    
  vector[I] beta_ij;     // vector of interactions which are inferrable by the NDDM
  
  vector<lower=0>[1] response1; // first species-specific response parameter
  // constrained to be positive for identifiability
  vector[S-1] responseSm1; // other species-specific response parameters

  unit_vector[K] effect; // species-specific effect parameter

} 

transformed parameters {
  
  // transformed parameters constructed from parameters above
  vector[S] response;        // combined vector of species-specific responses
  
  vector[N] mu;              // the RIM linear predictor for perform (here seed production)
  vector[N] mu2;             // the NDDM linear predictor for perform (here seed production)
  
  matrix[S, K] ri_betaij;   // interaction matrix for the RIM
  matrix[S, K] ndd_betaij;  // interaction matrix for the NDDM
  
  
  // stitch together the response values
  response = append_row(response1, responseSm1);
  // get RIM interactions
  ri_betaij = response*effect';
   
  // RIM estimates all interactions
  for(n in 1:N) {
       mu[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], ri_betaij[species_ID[n], ]));  
  }
  
  // NDDM estimates inferrable interactions, and uses RIM estimates when non-inferrable:
  
  ndd_betaij = ri_betaij; // initialise nddm interaction matrix to rim estimates
  // match inferrable interactions parameters to the correct position in the interaction matrix
  for(s in 1:S) {
    for(i in istart[s]:iend[s]) {
      ndd_betaij[irow[i], icol[i]] = beta_ij[i];
   }
  }
  // estimate inferrable interactions
  for(n in 1:N) {
        mu2[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], ndd_betaij[species_ID[n], ]));  
   }
   
} 

model {

  // priors
  beta_i0 ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
  beta_ij ~ normal(0,1);    // prior for interactions inferred by the NDDM
  
  response1 ~ normal(0, 1);   // constrained by parameter defintion to be positive
  responseSm1 ~ normal(0,1);
  // no prior needed for effect as we can use the default prior for the unit_vector

  // maximise likelihood for the RIM 
  for(n in 1:N) {
    perform[n] ~ neg_binomial_2(mu[n], (disp_dev[species_ID[n]]^2)^(-1));
    // in our case study, the response variable (seed production) shows a better fit to 
    // negative binomial (modify as necessary)
  }

  // add likelihood for the NDDM over inferrable interactions only.
  for (n in 1:N) {
    target += neg_binomial_2_lpmf(perform[n] | mu2[n], (disp_dev[species_ID[n]]^2)^(-1));
  }
  
}


