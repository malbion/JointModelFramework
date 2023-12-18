# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

library(rethinking)
library(reshape2)

# load required functions
source('data_prep.R')
source('simul_data.R')

# load simulated data 
set.seed(54)
simdat <- simul_data(S=3, T=3, p=0.4)   # p = proportion of interactions which are NOT observed 
df <- simdat[[1]]
sim_gamma <- simdat[[2]]
sim_interactions <- simdat[[3]]
# NB: if using real data or named species, ensure they are ordered alphabetically in the dataset

# identify focal and neighbouring species to be matched to parameter estimates
focalID <- unique(df$focal)  # this should return the names of unique focal groups in the order
# in which they are encountered in the dataframe - must be alphabetical
neighbourID <- colnames(df[ , -c(1:2)]) # should be ordered focal first (alphabetically), then
# non-focals in alphabetical order

# ensure neighbours are linearly independent across the whole dataset (see S1.2)
N_all <- df[ , neighbourID]
N_all <- apply(N_all, c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(df$focal)), N_all)
R_all <- pracma::rref(X_all)
Z_all <- t(R_all) %*% R_all
indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
  ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
}) #
all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep == 1)) warning('WARNING neighbours are not linearly independent') 


# prepare the data into the format required by STAN and the model code
stan.data <- data_prep(perform = 'seeds', 
                       focal = 'focal', 
                       nonNcols = 2, # number of columns that aren't neighbour abundances
                       df = df)


message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
message(paste0('Number of focal groups = ', length(focalID)))
message(paste0('Number of neighbour groups = ', length(neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$`T`)))

# Run the model! 
stan.seed <- 1234
fit <- stan(file = 'joint_model.stan', 
            data =  stan.data,               # named list of data
            chains = 4,
            warmup = 2000,          # number of warmup iterations per chain
            iter = 3000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10,
                           adapt_delta = 0.9), 
            seed = stan.seed
)

# check convergence
print(summary(fit, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit, pars=c("gamma_i","ndd_betaij"))
rstan::stan_rhat(fit, pars=c("gamma_i","ndd_betaij"))

# Get the full posteriors 
joint.post.draws <- extract.samples(fit)

# Select parameters of interest
param.vec <- fit@model_pars[!fit@model_pars %in% c('lp__')]

# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()
p.samples <- sapply(param.vec[!param.vec %in% c('ri_betaij', 'ndd_betaij')], function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 100)
  })  # this only works for parameters which are vectors
})


# WARNING: in the STAN model for annual wildflowers, parameter 'gamma_i' lies within an exponential,
# 'gamma_i' estimates must thus be exponentiated to return estimates of intrinsic performance
intrinsic.perf <- exp(p.samples$gamma_i)
colnames(intrinsic.perf) <- focalID

# Get interaction estimates
#--------------------------
inter_mat <- aperm(joint.post.draws$ndd_betaij, c(2, 3, 1))
rownames(inter_mat) <- focalID
colnames(inter_mat) <- neighbourID

# inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
# inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
# apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
# mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
# summary(fit), which is calculated from the full posterior distribution) 

# Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
# if a population dynamics model is used) in order to return per capita interaction strengths. 

