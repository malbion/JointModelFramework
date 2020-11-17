# Function to create a simulated dataset 


S <- 20   # number of focal groups / species
K <- 20   # number of neighbour focals 
pI <- 0.25   # proportion of interactions which are NOT observed 


# Neighbourhood abundance counts
#-------------------------------
# we assume we have a different number of observations for each focal
S_obs <- c(rpois(.15*S, 100), rpois(.15*S, 80), rpois(.5*S, 50), rpois(.2*S, 35))
K_Nmat <- matrix(data = NA, nrow = sum(S_obs), ncol = K)


# we assume that S = K, and that the number of observations of K is relative to their abundance
approx_K_obs <- round(jitter(10*(S_obs), amount = 30))
summedSK <- round(jitter(sapply(S_obs, '*', approx_K_obs)/1000, 20))

# randomly select observations that will be set to 0
obs_to_rm <- sample(seq_along(1:(S*K)), pI*S*K)

for(i in 1:(pI*S*K)) {
  summedSK[[obs_to_rm[i]]] <- 0
}

summedSK[1, ]
S_obs[1]



# Parameters
#-----------

# species-specific intrinsic performance 
sim_a <- rnbinom(20, mu = 70, size = 1)  

# 'true' interaction strengths 
a <- rnorm( 8000, 0, .3 )
b <- rnorm( 1000, -0.3, .2 )
c <- rnorm( 2000,  0.7, .3 )
customdist <- c( a, b, c ) ## make a weird dist with Kurtosis and Skew
hist(customdist, freq=FALSE)
mean(customdist)
sim_truealpha <- matrix(data = sample(customdist, S*K, replace = F),
                        nrow = S, ncol = K)



