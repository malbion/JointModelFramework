# This function takes a dataframe and returns a list (stan.data)
# with all the information necessary for STAN to run the model

# The dataframe should be organised with 1 observation of individual performance per row 
# Columns should include the measure of performance (perform), a column identifying which 
# species or element that observation belongs to (focal), and a column for each
# potential neighbouring species or element counting their abundance for that observation
# e.g.: 

#   focal  |  seeds  | species1 | species2 | species3
# ---------------------------------------------------
# species1 |   5     |    0     |    2     |    1
# species1 |   2     |    3     |    1     |    0
# species2 |   10    |    1     |    1     |    4

  
  
data_prep <- function(perform = "seeds", # column name for performance indicator
                      focal = 'focal',   # column name for species or element indicator
                      nonNcols = 2,      # number of columns before neighbour abundances begin
                      df = NULL,         # dataframe object
                      ...){
  

  # set up the data in list format as preferred by STAN:
  stan.data <- list()
  
  # create a matrix which tallies the number of realised interactions for each focal and neighbour
  # here we define
  counts <- df[ , -c(1:nonNcols)] # keep only neighbour abundances
  counts[counts>0] <- 1 # here we de
  counts <- split(counts, as.factor(df[ , focal]))
  obs <- do.call(rbind, lapply(counts, colSums))
  
  # INTEGERS
  stan.data$S <- length(unique(df[ , focal]))  # number of focal elements (species)
  stan.data$N <- nrow(df)                      # number of observations
  stan.data$K <- ncol(df[ , -c(1:nonNcols)])   # number of neighbours
  stan.data$I <- length(obs[obs>0])            # number of realised interactions

  
  # VECTORS
  stan.data$species_ID <- as.numeric(as.factor(df[ , focal])) 
  stan.data$perform <- df[ , perform]
  
  # set up indices to place realised interactions in the interaction matrix
  # first count the number of interactions realised for each focal species
  stan.data$inter_per_species <- obs
  stan.data$inter_per_species[stan.data$inter_per_species > 0] <- 1 # this counts every interaction
  # for which a focal i and neighbour j cooccur at least once as realised. 
  stan.data$inter_per_species <- rowSums(stan.data$inter_per_species)
  # column index in the interactions matrix for each realised interaction
  stan.data$icol <- unlist(apply(ifelse(obs > 0, T, F), 1, which))
  names(stan.data$icol) <- NULL
  stan.data$icol <- as.vector(stan.data$icol)
  # begin the row index
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]])
  # begin the start and end indices for the vector of realised interactions per species
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
  
  # MODEL MATRIX
  stan.data$X <- as.matrix(df[ , -c(1:nonNcols)]) 
  
  # Done!
  return(stan.data)
  
}