# This function takes a dataframe and returns a list (stan.data)
# with all the information necessary for STAN to run the model

# The dataframe should be organised with 1 observation of individual performance per row 
# Columns should include the measure of performance (perform), a column identifying which 
# species or component/group that observation belongs to (focal), and a column for each
# potential neighbouring species or component counting their abundance for that observation
# e.g.: 

#   focal  |  seeds  | species1 | species2 | species3
# ---------------------------------------------------
# species1 |   5     |    0     |    2     |    1
# species1 |   2     |    3     |    1     |    0
# species2 |   10    |    1     |    1     |    4

  
  
data_prep <- function(perform = "seeds", # column name for performance indicator
                      focal = 'focal',   # column name for species or component indicator
                      nonNcols = 2,      # number of columns before neighbour abundances begin
                      df = NULL,       # dataframe object
                      ...){
  

  # set up the data in list format as preferred by STAN:
  stan.data <- list()
  
  # create a matrix which tallies the number of observations for each focal and neighbour
  counts <- df[ , -c(1:nonNcols)] # keep only neighbour abundances
  counts[counts>0] <- 1
  counts <- split(counts, as.factor(df[ , focal]))
  obs <- do.call(rbind, lapply(counts, colSums))
  
  # integers
  stan.data$S <- length(unique(df[ , focal]))  # number of focal components (species)
  stan.data$N <- nrow(df)                      # number of observations
  stan.data$K <- ncol(df[ , -c(1:nonNcols)])   # number of neighbours
  stan.data$I <- length(obs[obs>0])              # number of observed interactions

  
  # vectors
  stan.data$species_ID <- as.numeric(as.factor(df[ , focal]))
  stan.data$perform <- df[ , perform]
  
  
  # set up indices to place observed interactions in the alpha matrix
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
  stan.data$X <- as.matrix(df[ , -c(1:nonNcols)]) 
  
  # Number of observations per interaction - for model weighting?
  stan.data$Obs <- as.vector(apply(obs, 1, c))   # vector of the number of observations for each interactions
  stan.data$Obs <- stan.data$Obs[stan.data$Obs>0] # remove unobserved interactions
  
  # Done!
  return(stan.data)
  
}