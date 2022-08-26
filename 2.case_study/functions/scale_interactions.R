# Scale the interaction parameters into betas
# Returns a 3D array of dimensions [N_species, N_neighbs, N_samples]

scale_interactions <- function(betas,
                               lambdas,   # growth.rates.samples
                               key_speciesID,
                               key_neighbourID,
                               comm,
                               ...){ 
  
  N_samples <- dim(betas)[1]
  # order as focals / neighbours / samples
  betamat <- aperm(betaij, c(2, 3, 1))
  
  # average seed rates 
  seeds <- read.csv('data/seed_rates.csv', row.names = 1)
  # Scale lambdas into growth rates including the seed rates (demographic param)
  # r_i = ln(eta) = ln( (g_i * lambda_i ) / (1 - (1 - g_i) s_i) )
    r_i <- sapply(c(1:length(key_speciesID)), function(x) {
    prod <- seeds[key_speciesID[x], 'germ'] / (1 - (1 - seeds[key_speciesID[x], 'germ'])*seeds[key_speciesID[x], 'surv'])
    if (prod < 0.1) {prod <- 0.1}
      log(lambdas[,x ] * prod)
    })
    
    dimnames(r_i) <- list('samples' = seq(1:N_samples), 'species' = key_speciesID)
    ri_means <- colMeans(r_i)[order(colMeans(r_i), decreasing = T)]
    message(paste0('Growth rates for community ', comm))
    print(ri_means)
    
    #----------------- Removing growth rates < 1 ?? ------------------------------------------#
    if(length(ri_means[ri_means < 1]) > 0){
      # If any growth rates are below 1 .... 
      to.rplc <- names(ri_means[ri_means < 1])
      # replace with the next lowest growth rate (that's over 1)
      rplc.ri <- min(ri_means[ri_means > 1])
      rplc.sp <- names(ri_means[ri_means == rplc.ri])
      
      message(paste0('Growth rates for ', to.rplc, ' replaced by ', rplc.sp, '\n'))
      
      for(i in to.rplc) {
        r_i[ , i] <- r_i[, rplc.sp]
      }
 }
  
    #-------------------------------------------------------------------------------------#
    write.csv(r_i, paste0('model/transformed/ri_demogrwth.csv'), row.names = T)
    
    
    # Scale the betas by dividing by growth rate
    # scaled beta = model_param/r_i = model_param/ln(eta)
    # reminder: model_param == alpha*g (in plant pop dyn model)
    beta_mat <- sapply(c(1:length(key_speciesID)), function(x) {
      # sample each interaction from the community matrix
      bij_samples <- t(betamat[x, , ])
      # sample the growth rates
      ri_samples <- r_i[ , x]
      # divide each row - aij_samples must be [N_samples, N_neighbours]
      bij_samples/ri_samples
    }, simplify = 'array', USE.NAMES = T)
    
    # SAMPLE from the 80% CI for each tranformed beta
    beta_df <- as.data.frame(beta_mat)
    beta_CI <- apply(beta_df, 2, function(param) {
      sample(param[param > quantile(param, 0.1) & 
                     param < quantile(param, 0.9)], size = 1000)
    })
    # get betas back into a focals / neighbours / samples format
    beta_scaled <- array(data = unlist(beta_CI), 
                          dim = c(1000, length(key_neighbourID), length(key_speciesID)), 
                          dimnames = list('samples' = seq(1, 1000), 'neighbour' = key_neighbourID, 'species' = key_speciesID))
    beta_scaled <- aperm(beta_scaled, c(3,2,1))
    
    dimnames(beta_scaled)[[1]] <- key_speciesID
    names(dimnames(beta_scaled))[1] <- 'species'
    dimnames(beta_scaled)[[3]] <- seq(1:1000)
    names(dimnames(beta_scaled))[3] <- 'samples'
    attributes(dimnames(beta_scaled)[[1]]) <- NULL    
    attributes(dimnames(beta_scaled)[[2]]) <- NULL

  # Check that the rows and columns are properly ordered
  if(identical(dimnames(beta_scaled)[[1]], dimnames(beta_scaled)[[2]][1:dim(beta_scaled)[[1]]]) == F) {
    message('WARNING: rownames and column names do no match')
  }
  
  return(beta_scaled)
  
  }


