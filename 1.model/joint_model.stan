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
  
  vector[S] a;    // species-specific intercept, i.e. log(lambda) (intrinsic fitness)
  vector<lower=0>[S] disp_dev; // species-specific param for dispersion deviation, 
  // disp_dev = 1/sqrt(phi) 
  
  
  vector[I] interactions;     // vector of interactions which have been observed
  vector[K] effect;            // competitive effect of neighbours, same across all focals
  // can be facilitative (-) or competitive (+)
  vector<lower=0>[S] response; // species-specific competitive response parameter
  // >= 0 to avoid bimodality in response and effect  

  real<lower=0> sigma_alph;	// variance for the ifm alphas
} 

transformed parameters {
  
  // transformed parameters constructed from params above
  vector[N] mu;              // the linear predictor
  matrix[S, K] ifm_alpha; 	 // community matrix 
  vector[I] re;              // interactions as calculated by the re model
  
  ifm_alpha = rep_matrix(0, S, K);   // fill the community matrix with 0 (instead of NA)
  
  // match observed interactions parameters to the correct position in the community matrix
  for(s in 1:S) {
    for(i in istart[s]:iend[s]) {
      ifm_alpha[irow[i], icol[i]] = interactions[i];
    }
  }
  
  // individual fitness model 
  for(n in 1:N) {
       mu[n] = exp(a[species_ID[n]] - dot_product(X[n], ifm_alpha[species_ID[n], ]));  
  }
  
  // build a vector of interaction parameters based on the response effect model approach
  for (i in 1:I) {
    re[i] = response[irow[i]]*effect[icol[i]];
  }
  
} 

model {

  // priors
  a ~ cauchy(0,10); // prior for the intercept following Gelman 2008
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
  response ~ normal(0, 1);   // 
  effect ~ normal(0, 1);
  sigma_alph ~ cauchy(0, 1);

  // seed production 
  for(n in 1:N) {
    seeds[n] ~ neg_binomial_2(mu[n], (disp_dev[species_ID[n]]^2)^(-1));
  }

  // response-effect interactions
  for (i in 1:I) {
    target += logistic_lpdf(re[i] | interactions[i], sigma_alph);
  }
  
}


