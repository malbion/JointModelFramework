/* NDDM & RIM joint model

code to run on all focal species at once

Bimler 2021
*/ 
  
  
data {
  int<lower=1> S;          // number of species (elements) 
  int<lower=1> N;          // number of observations (rows in model matrix)
  int<lower=0> K;          // number of neighbours (columns in model matrix)
  int<lower=0> I;          // number of realised interactions 
  
  int<lower=0> species_ID[N];   // index matching species to observations
  int<lower=0> perform[N];      // response variable 
  
  int<lower=0> istart[S];       // indices matching species to realised interactions
  int<lower=0> iend[S];
  int<lower=0> icol[I];
  int<lower=0> irow[I];
  
  matrix[N,K] X;         // neighbour abundances (the model matrix)
  matrix[S, K] Q;	  // inferrable interactions

} 

parameters {
  
  vector[S] beta_i0;    // species-specific intercept 
  vector<lower=0>[S] disp_dev; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)

  // vector[I] beta_ij;     // vector of interactions which have been realised
  
  vector<lower=0>[1] response1; // first species-specific response parameter
  // constrained to be positive for identifiability
  vector[S-1] responseSm1; // other species-specific response parameters

  unit_vector[K] effect; // species-specific effect parameter

  real<lower=0> sigma;	// scale parameter for the logistic distribution (used to estimate re's)
} 

transformed parameters {
  
  // transformed parameters constructed from parameters above
  vector[S] response;        // combined vector of species-specific responses
  // matrix[S, K] inter_mat;    // the community interaction matrix
  // vector[I] re;              // interactions as calculated by the re model
  vector[N] mu;              // the linear predictor for perform (here seed production)
  matrix[S, K] rim_betaij;
  
  // inter_mat = rep_matrix(0, S, K); // fill the community interaction matrix with 0 (instead of NA)
  
  // stitch together the response values
  response = append_row(response1, responseSm1);
  // get RIM interactions
  rim_betaij = response*effect';
  
  // response-impact model runs on all interactions
  for(n in 1:N) {
       mu[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], rim_betaij[species_ID[n], ]));  
  }
  
  
  // match observed interactions parameters to the correct position in the community matrix
  // for(s in 1:S) {
  //   for(i in istart[s]:iend[s]) {
  //     inter_mat[irow[i], icol[i]] = beta_ij[i];
  //  }
  // }
  
  // neighbour density-dependent model 
  // for(n in 1:N) {
  //      mu[n] = exp(beta_i0[species_ID[n]] - dot_product(X[n], inter_mat[species_ID[n], ]));  
  // }


  
  // build a vector of interaction parameters based on the response impact model
  // for (i in 1:I) {
  //   re[i] = response[irow[i]]*effect[icol[i]];
  // }
} 

model {

  // priors
  beta_i0 ~ cauchy(0,10);   // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
  response1 ~ normal(0, 1);   // constrained by parameter defintion to be positive
  responseSm1 ~ normal(0,1);
  sigma ~ cauchy(0, 1);
  // no prior needed for effect as we can use the default prior for the unit_vector

  // seed production, i.e. P_i
  for(n in 1:N) {
    perform[n] ~ neg_binomial_2(mu[n], (disp_dev[species_ID[n]]^2)^(-1));
    // in our case study, seed production shows a better fit to a negative binomial 
    // than poisson distribution
  }

  // response-effect interactions
  // for (i in 1:I) {
  //   target += logistic_lpdf(re[i] | beta_ij[i], sigma);
  // }
  
}


