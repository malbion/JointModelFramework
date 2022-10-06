/* NDDM 

code to run on all focal species at once

Bimler 2022
*/ 
  
  
data {
  int<lower=1> S;          // number of species (elements) 
  int<lower=1> N;          // number of observations (rows in model matrix)
  int<lower=0> T;          // number of neighbours (columns in model matrix)
  int<lower=0> I;          // number of identifiable interactions 
  
  int<lower=0> species_ID[N];   // index matching species to observations
  int<lower=0> perform[N];      // response variable 
    
  int<lower=0> istart[S];       // indices matching species to identifiable interactions
  int<lower=0> iend[S];
  int<lower=0> icol[I];
  int<lower=0> irow[I];
  
  matrix[N,T] X;         // neighbour abundances (the model matrix)

} 

parameters {
  
  vector[S] gamma_i;    // species-specific intercept 
    
  vector<lower=0>[S] disp_dev; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
    
  vector[I] beta_ij;     // vector of interactions which can be identified
  
} 

transformed parameters {
  
  // transformed parameters constructed from parameters above

  vector[N] mu;              // the linear predictor for perform (here seed production)
  matrix[S, T] ndd_betaij;  // interaction matrix for the NDDM

  ndd_betaij = rep_matrix(0, S, T); // fill the community interaction matrix with 0 (instead of NA)
    
  // match observed interactions parameters to the correct position in the community matrix
  for(s in 1:S) {
    for(i in istart[s]:iend[s]) {
      ndd_betaij[irow[i], icol[i]] = beta_ij[i];
   }
  }
  // estimate identifiable interactions
  for(n in 1:N) {
        mu[n] = exp(gamma_i[species_ID[n]] - dot_product(X[n], ndd_betaij[species_ID[n], ]));  
   }

 
} 

model {

  // priors
  gamma_i ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
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



