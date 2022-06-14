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
  
  # MATRIX OF INFERRABLE INTERACTIONS
  # this is done species by species
  Q <- t(sapply(levels(as.factor(df[ , focal])), function(f){
    
    N_i <- as.matrix(df[df[ , focal] == f, -c(1:nonNcols)])
    X_i <- cbind(1,N_i)
    R_i <- pracma::rref(X_i)
    Z_i <- t(R_i) %*% R_i
    
    # param k is inferrable if its corresponding row/column is all 0 except for the k'th element
    # ignore intercept because we always want to include it
    sapply(seq(2, dim(Z_i)[1], 1), function(k){ 
      ifelse(Z_i[k, k] == 1 & sum(Z_i[k, -k]) == 0, 1, 0)
    }) 
    
  }))
  stan.data$Q <- Q
  # Q is a matrix of focal x neighbours, if Q[i, j] = 1 then the interaction between i and j is inferrable
  
  # INTEGERS
  stan.data$S <- length(unique(df[ , focal]))  # number of focal elements (species)
  stan.data$N <- nrow(df)                      # number of observations
  stan.data$K <- ncol(df[ , -c(1:nonNcols)])   # number of neighbours
  stan.data$I <- sum(Q)                        # number of inferred interactions

  
  # VECTORS
  stan.data$species_ID <- as.numeric(as.factor(df[ , focal])) 
  stan.data$perform <- df[ , perform]
  
  # set up indices to place realised interactions in the interaction matrix
  # first count the number of identifiable interactions for each focal species
  stan.data$inter_per_species <- rowSums(Q)
  # column index in the interactions matrix for each inferrable interaction
  stan.data$icol <- unlist(apply(ifelse(Q > 0, T, F), 1, which))
  names(stan.data$icol) <- NULL
  stan.data$icol <- as.vector(stan.data$icol)
  # begin the row index
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]])
  # begin the start and end indices for the vector of inferrable interactions per species
  stan.data$istart <- 1
  stan.data$iend <- stan.data$inter_per_species[[1]] 
  
  # populate indices for all the other species
  for(s in 2:stan.data$S) {
    # starting position of a_ij's for i in the vector of inferrable interactions (ie the 1st 'j')
    stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
    # end position of a_ij's for i in the vector of inferrable interactions (the last 'j')
    stan.data$iend[s] <-  sum(stan.data$inter_per_species[1:s])
    # row index in the alpha matrix for each inferrable interaction
    stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
  }
  
  # MODEL MATRIX
  stan.data$X <- as.matrix(df[ , -c(1:nonNcols)]) # neighbour abundances
  
  
  # Done!
  return(stan.data)
  
}