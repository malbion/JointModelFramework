# Function to create a simulated dataset 

simul_data <- function(
                        # number of species has to be a multiple of 10
                        S = 10,   # number of focal groups / species
                        K = 10,   # number of neighbour focals 
                        pI = 0.25,   # proportion of interactions which are NOT observed 
                        ...
) {
  
  # Neighbourhood abundance counts
  #-------------------------------
  # we assume we have a different number of observations for each focal
  S_obs <- c(rpois(.1*S, 100), rpois(.2*S, 80), rpois(.5*S, 50), rpois(.2*S, 35))
  
  # we assume that S = K, and that the number of observations of K is relative to their abundance
  approx_K_obs <- round(jitter(10*(S_obs), amount = 30))
  summedSK <- round(jitter(sapply(S_obs, '*', approx_K_obs)/1000, 20))
  
  # randomly select observations that will be set to 0
  obs_to_rm <- sample(seq_along(1:(S*K)), pI*S*K)
  # remove unobserved interactions
  for(i in 1:(pI*S*K)) {
    summedSK[[obs_to_rm[i]]] <- 0
  }
  
  # set up an empty matrix with columns for each neighbours
  K_Nmat <- matrix(ncol = K)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    foo <- matrix(data = 0, nrow = S_obs[s], ncol = K)
    for(j in 1:K) {
      if(summedSK[s, j] > 0)
        for (i in 1:summedSK[s, j]){ 
          # randomly select an observation
          cellN <- round(runif(1, 1, S_obs[s]))
          # fill it with a 1 
          foo[cellN, j] <- foo[cellN, j] + 1  
          # and so on until neighbours are all accounted for
        }
    }
    if(identical(summedSK[s , ], colSums(foo)) == F) message('Error in abundance tallies!')
    K_Nmat <- rbind(K_Nmat, foo)
  }
  K_Nmat <- K_Nmat[-1, ]  # remove the NA row
  if (nrow(K_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  if (sum(colSums(K_Nmat) != colSums(summedSK)) > 1) message('Error in neighbour abundances')
  colnames(K_Nmat) <- paste0('K', 1:K)
  
  # Parameters
  #-----------
  # log of species-specific intrinsic performance
  sim_a <- runif(S, 2, 7)  
  
  # 'true' interaction strengths 
  a <- rnorm( 8000, 0, .3 )
  b <- rnorm( 2000, -0.7, .3 )
  c <- rnorm( 1000,  0.3, .2 )
  customdist <- c( a, b, c ) ## make a weird dist with Kurtosis and Skew
  sim_truealpha <- matrix(data = sample(customdist, S*K, replace = F),
                          nrow = S, ncol = K)
  
  # Observed seed set
  #------------------
  seeds <- rep(NA, length = sum(S_obs))
  counter <- 0
  for(s in 1:S) {
    for (i in (counter+1):(counter + S_obs[s])) {
      # multiply neighbour abundances by 'true' alpha values
      seeds[i] <- round(exp(sim_a[s] - sum(K_Nmat[i, ] * sim_truealpha[s, ])))
    }
    counter <- counter + S_obs[s]
  } 
  
  # Simulated dataset
  #------------------
  focal <- do.call('c', mapply(rep, 1:S, S_obs, SIMPLIFY = F))    # focal identifier
  simdata <- cbind(focal, seeds, K_Nmat)
  simdata <- as.data.frame(simdata)
  
  return(list(simdata, sim_a, sim_truealpha))
  
}




  
