# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

library(rethinking)
library(reshape2)
# library(coda)

# load required functions
source('data_prep.R')
source('return_inter_array.R')
source('simul_data.R')

# load simulated data 
simdat <- simul_data(S=10, K=10, p=0.25)
df <- simdat[[1]]
sim_a <- simdat[[2]]
sim_interactions <- simdat[[3]]

# identify focal and neighbouring species to be matched to parameter estimates
focalID <- unique(df$focal)  # this should return the names of unique focal groups in the order
# in which they are encountered in the dataframe
neighbourID <- colnames(df[ , -c(1:2)])

# ensure neighbours are linearly independent across the whole dataset
N_all <- df[ , neighbourID]
N_all <- apply(N_all, c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(df$focal)), N_all)
R_all <- pracma::rref(X_all)
Z_all <- t(R_all) %*% R_all
indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
  ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
}) #
all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep == 1)) message('WARNING neighbours are not linearly independent') 


# prepare the data into the format required by STAN and the model code
stan.data <- data_prep(perform = 'seeds', 
                       focal = 'focal', 
                       nonNcols = 2, # number of columns that aren't neighbour abundances
                       df = df)


message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
message(paste0('Number of focal groups = ', length(focalID)))
message(paste0('Number of neighbour groups = ', length(neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$K)))

# Run the model! 
fit <- stan(file = 'joint_model.stan',
            data =  stan.data,               # named list of data
            chains = 1,
            warmup = 100,          # number of warmup iterations per chain
            iter = 600,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10,
                           adapt_delta = 0.95)
)

# As well as the usual traceplots etc., convergence one 1 chain can be checked 
# using the geweke.diag() function from the coda package, e.g.:
matrix_of_draws <- as.matrix(fit)
gew <- geweke.diag(matrix_of_draws)

# Get the full posteriors 
joint.post.draws <- extract.samples(fit)

# Select parameters of interest
param.vec <- c('beta_i0', 'beta_ij', 'effect', 'response', 're', 'inter_mat',
               'mu', 'disp_dev', 'sigma') 


######### diagnostics only ####################
source('../2.case_study/functions/stan_modelcheck_rem.R')
detach('package:coda')  # maybe this isnt necessary after not loading it directly anymore

stan_diagnostic(fit, 'output')
stan_model_check(fit, 'output', params = param.vec)
stan_post_pred_check(joint.post.draws, 'output', stan.data)

log_post <- unlist(extract(fit, 'lp__'))
write.csv(log_post, file = paste0('output/log_post.csv'), row.names = F)

##################################################################################

# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()
p.samples <- sapply(param.vec[param.vec != 'sigma' & param.vec != 'inter_mat'], function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  })  # this only works for parameters which are vectors
})
# there is only one sigma_alph parameter so we must sample differently:
p.samples[['sigma']] <- sample(joint.post.draws$sigma[
  joint.post.draws$sigma > quantile(joint.post.draws$sigma, 0.1) & 
    joint.post.draws$sigma < quantile(joint.post.draws$sigma, 0.9)], size = 1000)

# WARNING: in the STAN model for annual wildflowers, parameter 'a' lies within an exponential,
# 'a' estimates must thus be exponentiated to return estimates of intrinsic performance
intrinsic.perf <- exp(p.samples$beta_i0)
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

