# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

library(rethinking)
library(reshape2)

# load required functions
source('data_prep.R')
source('return_inter_array.R')
source('simul_data.R')

# load the data 
simdat <- simul_data(S=10, K=10, p=0.25)
df <- simdat[[1]]
sim_a <- simdat[[2]]
sim_interactions <- simdat[[3]]

# prepare the data into the format required by STAN and the model code
stan.data <- data_prep(perform = 'seeds', 
                       focal = 'focal', 
                       nonNcols = 2, # number of columns that aren't neighbour abundances
                       df = df)

# identify focal and neighbouring species to be matched to parameter estimates
focalID <- unique(df$focal)  # this should return the names of unique focal groups in the order
# in which they are encountered in the dataframe
neighbourID <- colnames(df[ , -c(1:2)])

message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
message(paste0('Number of focal groups = ', length(focalID)))
message(paste0('Number of neighbour groups = ', length(neighbourID)))


# Run the model! 
fit <- stan(file = 'joint_model.stan',
            data =  stan.data,               # named list of data
            chains = 4,
            warmup = 1000,          # number of warmup iterations per chain
            iter = 3000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10,
                           adapt_delta = 0.9)
)

# Get the full posteriors 
joint.post.draws <- extract.samples(fit)

# Select parameters of interest
param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat',
               'mu', 'disp_dev', 'sigma_alph')

# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()
p.samples <- sapply(param.vec[param.vec != 'sigma_alph' & param.vec != 'inter_mat'], function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  })  # this only works for parameters which are vectors
})
# there is only one sigma_alph parameter so we must sample differently:
p.samples[['sigma_alph']] <- sample(joint.post.draws$sigma_alph[
  joint.post.draws$sigma_alph > quantile(joint.post.draws$sigma_alph, 0.1) & 
    joint.post.draws$sigma_alph < quantile(joint.post.draws$sigma_alph, 0.9)], size = 1000)

# WARNING: in the STAN model for annual wildflowers, parameter 'a' lies within an exponential,
# 'a' estimates must thus be exponentiated to return estimates of intrinsic performance
intrinsic.perf <- exp(p.samples$a)
colnames(intrinsic.perf) <- focalID

# Build matrices of interaction estimates
#----------------------------------------
inter_mat <- return_inter_array(joint.post.draws, 
                                response = p.samples$response,
                                effect = p.samples$effect,
                                focalID,
                                neighbourID)
# inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
# inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
# apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
# mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
# summary(fit), which is calculated from the full posterior distribution) 

# Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
# if a population dynamics model is used) in order to return per capita interaction strengths. 

