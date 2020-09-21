/* 2018_CompNet
model

code to run a joint IF and RE model on 
all focal species at once
*/ 
  
  
  data {
    int<lower=1> S;          // number of species (rows in model matrix)
    int<lower=1> N;          // number of observations
    int<lower=0> K;          // number of neighbours (columns in model matrix)
    int<lower=0> I;          // number of interactions observed
    
    int<lower=0> species_ID[N];   // species index for observations
    int<lower=0> seeds[N];        // response variable 
    
    int<lower=0> istart[S];       // indices matching species to interactions
    int<lower=0> iend[S];
    int<lower=0> icol[I];
    int<lower=0> irow[I];
    
    matrix[N,K] X;         // neighbour abundances, the model matrix
  
  } 

parameters {
  
  vector[S] a;    // species-specific intercept, i.e. log(beta_i_0) 
  vector<lower=0>[S] disp_dev; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (seeds)
  // disp_dev = 1/sqrt(phi) 
  
  
  vector[I] beta_ij;     // vector of interactions which have been observed
  vector[K] effect;            // competitive effect of neighbours, same across all 
  // focals, can be facilitative (-) or competitive (+)
  vector<lower=0>[S] response; // species-specific competitive response parameter
  // >= 0 to avoid bimodality in response and effect  

  real<lower=0> sigma;	// scale parameter for the logistic distribution (used to estimate re's)
} 

transformed parameters {
  
  // transformed parameters constructed from paramaters above
  vector[N] mu;              // the linear predictor for seed production 
  matrix[S, K] inter_mat;     // the community interaction matrix 
  vector[I] re;              // interactions as calculated by the re model
  
  inter_mat = rep_matrix(0, S, K); // fill the community interaction matrix with 0 (instead of NA)
  
  // match observed interactions parameters to the correct position in the community matrix
  for(s in 1:S) {
    for(i in istart[s]:iend[s]) {
      inter_mat[irow[i], icol[i]] = beta_ij[i];
    }
  }
  
  // individual fitness model 
  for(n in 1:N) {
       mu[n] = exp(a[species_ID[n]] - dot_product(X[n], inter_mat[species_ID[n], ]));  
  }
  
  // build a vector of interaction parameters based on the response effect model 
  for (i in 1:I) {
    re[i] = response[irow[i]]*effect[icol[i]];
  }
} 

model {

  // priors
  a ~ cauchy(0,10); // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
  response ~ normal(0, 1);   // defining a lower limit aboves truncates the normal distribution
  effect ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // seed production, i.e. F_i
  for(n in 1:N) {
    seeds[n] ~ neg_binomial_2(mu[n], (disp_dev[species_ID[n]]^2)^(-1));
    // in our case study, seed production shows a better fit to a negative binomial than poisson   
    // distribution
  }

  // response-effect interactions
  for (i in 1:I) {
    target += logistic_lpdf(re[i] | beta_ij[i], sigma);
  }
  
}


