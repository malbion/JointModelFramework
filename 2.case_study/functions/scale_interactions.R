# Malyon Bimler
# CompNet 2019 

# Scale the interaction parameters into betas
# Returns a 3D array of dimensions [N_species, N_neighbs, N_samples]

scale_interactions <- function(betas,
                               lambdas,   # growth.rates.samples
                               key_speciesID,
                               key_neighbourID,
                               comm,
                               ...){ 
  
  N_samples <- nrow(betas)

  # Recreate community matrices (focals x neighbours, with 3rd dim = samples taken from posterior draws)
  alpha_mat <- array(data = unlist(betas), 
                     dim = c(N_samples, length(key_neighbourID), length(key_speciesID)), 
                     dimnames = list('samples' = seq(1, N_samples), 'neighbour' = key_neighbourID, 'species' = key_speciesID))
  alpha_mat <- aperm(alpha_mat, c(3,2,1))

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
        r_i[ , i] <- sample(r_i[, rplc.sp], length(r_i[, rplc.sp]), replace = F)
      }
 }
  
    #-------------------------------------------------------------------------------------#
    write.csv(r_i, paste0('model/transformed/ri_demogrwth.csv'), row.names = T)
    
    # Check all species names match! 
    ifelse(identical(dimnames(alpha_mat)[[1]], key_speciesID) == F, 
           message("ERROR! Comm matrix row names don't match species key"),
           ifelse(identical(colnames(r_i), key_speciesID) == F, 
                  message("ERROR! Growth rate column names don't match species key"), 
                  'Key species ID names match'))
    
    
    # Scale the betas by dividing by growth rate
    # scaled alpha = model_param/r_i = model_param/ln(eta)
    # reminder: model_param == alpha*g
    alpha_mat <- sapply(c(1:length(key_speciesID)), function(x) {
      # sample each interaction from the community matrix
      aij_samples <- alpha_mat[x, , ]
      aij_samples <- apply(aij_samples, 1, sample, N_samples, replace = T)
      # sample the growth rates
      ri_samples <- sample(r_i[ , x], N_samples, replace = T)
      # check dimensions 
      if(dim(aij_samples)[1] != length(ri_samples)) {
        message('Error!: aij samples have incorrect dim')}  # error code 
      # divide each row - aij_samples must be [N_samples, N_neighbours]
      aij_samples/ri_samples
    }, simplify = 'array', USE.NAMES = T)
    
    alpha_mat <- aperm(alpha_mat, c(3, 2, 1))
    dimnames(alpha_mat)[[1]] <- key_speciesID
    names(dimnames(alpha_mat))[1] <- 'species'
    dimnames(alpha_mat)[[3]] <- seq(1:N_samples)
    names(dimnames(alpha_mat))[3] <- 'samples'
    attributes(dimnames(alpha_mat)[[1]]) <- NULL

  # Check that the rows and columns are properly ordered
  if(identical(dimnames(alpha_mat)[[1]], dimnames(alpha_mat)[[2]][1:dim(alpha_mat)[[1]]]) == F) {
    message('WARNING: rownames and column names do no match')
  }
  
  return(alpha_mat)
  
  }


