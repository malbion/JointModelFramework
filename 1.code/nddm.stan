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
    
  vector[I] beta_ij;     // vector of interactions which can be inferred
  
} 

transformed parameters {
  
  // transformed parameters constructed from parameters above

  vector[N] mu;              // the linear predictor for perform (here seed production)
  matrix[S, K] ndd_betaij;  // interaction matrix for the NDD model

  ndd_betaij = rep_matrix(0, S, K); // fill the community interaction matrix with 0 (instead of NA)
    
  // match observed interactions parameters to the correct position in the community matrix
  for(s in 1:S) {
    for(i in istart[s]:iend[s]) {
      ndd_betaij[irow[i], icol[i]] = beta_ij[i];
   }
  }
  
  for(n in 1:N) {
        mu[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], ndd_betaij[species_ID[n], ]));  
   }

 
} 

model {

  // priors
  beta_i0 ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
  beta_ij ~ normal(0,1);
  
 
  // seed production, i.e. P_i
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



