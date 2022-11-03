# Function to create a simulated dataset 

simul_data <- function(
                        S,   # number of focal groups / species
                        T,   # number of neighbour focals 
                        pI,   # proportion of interactions which are NOT observed 
                        ...
) {
    # Parameters
  #-----------
  # log of species-specific intrinsic performance
  sim_a <- runif(S, 2, 7)  
  
  # 'true' interaction strengths 
  sim_responses <- runif(S,-0.1,1)
  sim_effects <- runif(T)
  sim_truealpha <- matrix(sim_responses%*%t(sim_effects), S, T)

  # Neighbourhood abundance counts
  #-------------------------------
  # we assume we have a different number of observations for each focal
  S_obs <- c(rpois(.1*S, 100), rpois(.2*S, 80), rpois(.5*S, 50), rpois(.2*S, 35))
  if (length(S_obs) < S) {S_obs <- c(rpois(S-length(S_obs), 100), S_obs)} # workaround for species number < 10
  
  # we assume that S = T, and that the number of observations of T is relative to their abundance
  approx_T_obs <- round(jitter(10*(S_obs), amount = 30))
  summedST <- round(jitter(sapply(S_obs, '*', approx_T_obs)/1000, 20))
  
  # randomly select observations that will be set to 0
  obs_to_rm <- sample(seq_along(1:(S*T)), pI*S*T)
  # remove unobserved interactions
  for(i in seq_along(obs_to_rm)) {
    summedST[[obs_to_rm[i]]] <- 0
  }
  
  # set up an empty matrix with columns for each neighbours
  T_Nmat <- matrix(ncol = T)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    foo <- matrix(data = 0, nrow = S_obs[s], ncol = T)
    for(j in 1:T) {
      if(summedST[s, j] > 0)
        for (i in 1:summedST[s, j]){ 
          # randomly select an observation
          cellN <- round(runif(1, 1, S_obs[s]))
          # fill it with a 1 
          foo[cellN, j] <- foo[cellN, j] + 1  
          # and so on until neighbours are all accounted for
        }
    }
    if(identical(summedST[s , ], colSums(foo)) == F) message('Error in abundance tallies!')
    T_Nmat <- rbind(T_Nmat, foo)
  }
  T_Nmat <- T_Nmat[-1, ]  # remove the NA row
  if (nrow(T_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  if (sum(colSums(T_Nmat) != colSums(summedST)) > 1) message('Error in neighbour abundances')
  colnames(T_Nmat) <- paste0('T', 1:T)
  
  # Observed seed set
  #------------------
  seeds <- rep(NA, length = sum(S_obs))
  counter <- 0
  for(s in 1:S) {
    for (i in (counter+1):(counter + S_obs[s])) {
      # multiply neighbour abundances by 'true' alpha values
      seeds[i] <- rpois(1,exp(sim_a[s] - sum(T_Nmat[i, ] * sim_truealpha[s, ])))
    }
    counter <- counter + S_obs[s]
  } 
  
  # Simulated dataset
  #------------------
  focal <- do.call('c', mapply(rep, 1:S, S_obs, SIMPLIFY = F))    # focal identifier
  simdata <- cbind(focal, seeds, T_Nmat)
  simdata <- as.data.frame(simdata)
  
  return(list(simdata, sim_a, sim_truealpha))
  
}




  
