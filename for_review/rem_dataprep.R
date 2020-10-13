# This function takes the fecundities dataframe and formats it 
# for STAN

rem_dataprep <- function(fecundities,
                          ...) {
  
  
  # set up the data in list format as preferred by STAN:
  stan.data <- list()
  
  # create a matrix which tallies the number of observations for each focal and neighbour
  fec <- fecundities[ , -c(1:4)]
  fec[fec>0] <- 1
  fec <- split(fec, fecundities$focal)
  obs <- do.call(rbind, lapply(fec, colSums))
  
  # integers
  stan.data$S <- length(unique(fecundities$focal))  # number of species
  stan.data$N <- nrow(fecundities)                  # number of observations
  stan.data$K <- ncol(fecundities[ , -c(1:4)])      # number of neighbours
  stan.data$I <- length(obs[obs>0])        # number of observed interactions
  # stan.data$Z <- stan.data$S*stan.data$K   # total number of possible interactions         
  
  # vectors 
  stan.data$species_ID <- as.numeric(as.factor(fecundities$focal))
  stan.data$seeds <- fecundities$seeds
  
  # set up indices to place observed interaction in the alpha matrix
  # first count the number of interactions observed for each focal species
  stan.data$inter_per_species <- obs
  stan.data$inter_per_species[stan.data$inter_per_species > 0] <- 1
  stan.data$inter_per_species <- rowSums(stan.data$inter_per_species)
  # column index in the alpha matrix for each observed interaction
  stan.data$icol <- unlist(apply(ifelse(obs > 0, T, F), 1, which)) 
  names(stan.data$icol) <- NULL
  stan.data$icol <- as.vector(stan.data$icol)
  # begin the row index
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]])
  # begin the start and end indices for the vector of interactions per species 
  stan.data$istart <- 1
  stan.data$iend <- stan.data$inter_per_species[[1]] #???
  # populate indices for all the other species
  for(s in 2:stan.data$S) {
    # starting position of a_ij's for i in the vector of observed interactions (ie the 1st 'j')
    stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
    # end position of a_ij's for i in the vector of observed interactions (the last 'j')
    stan.data$iend[s] <-  sum(stan.data$inter_per_species[1:s])
    # row index in the alpha matrix for each observed interaction
    stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
  }

  # model matrix 
  stan.data$X <- as.matrix(fecundities[ , -c(1:4)])  
  
  # Number of observations per interaction - for model weighting?
  stan.data$Obs <- as.vector(apply(obs, 1, c))   # vector of the number of observations for each interactions
  stan.data$Obs <- stan.data$Obs[stan.data$Obs>0] # remove unobserved interactions
  
  
  return(stan.data)
}



