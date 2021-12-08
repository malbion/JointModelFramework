# Objective: running a joint model combining IF and RE to estimate all interactions at the same time. 
# 
# Run in Terminal: 
#   
#   R CMD BATCH '--args 0' model.R model/0_model_script.Rout
# 
# substituting '0' for the comm number (0, 1, 2, or 3) 
# 
# This calls up model.R, which runs joint_model.stan. 
# results got to model/output, /transformed and /validation


# One model to rule them all: 
# 1. Estimate observed interactions with IFM 
# 2. Estimate response and effect from observed interactions

# Then: 
# 3. Extract the intrinsic growth rates
# 4. Estimate missing interactions using response and effect 
# 5. Scale the interactions into alphas

# PRELUDE
#-------
Sys.time()

# Get arguments from bash script
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# take comm name as an argument from bash
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
if (length(args)>1) {
  stop("Model can only be run on 1 comm at a time.n", call.=FALSE)
}
comm <- args[1]

# set up R environment
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

library(rethinking)
library(reshape2)
library(here)

# setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/3.case_study/')
setwd(here('3.case_study/'))


source('functions/rem_dataprep.R')
source('functions/stan_modelcheck_rem.R')
source('functions/scale_interactions.R')

# Load, verify & prepare data
#----------------------------
fecundities <- read.csv(paste0('data/fecundities', comm, '.csv'), stringsAsFactors = F)

# Keep note of how focals and neighbours are indexed
key_speciesID <- unlist(read.csv(paste0('data/key_speciesID', comm, '.csv'), stringsAsFactors = F))
key_neighbourID <- unlist(read.csv(paste0('data/key_neighbourID', comm, '.csv'), stringsAsFactors = F))


# ensure neighbours are linearly independent across the whole dataset
N_all <- apply(fecundities[ , 5:dim(fecundities)[2]], c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(fecundities$focal)), N_all)
R_all <- pracma::rref(X_all)
Z_all <- t(R_all) %*% R_all
indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
  ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
}) #
all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep == 1)) message('WARNING neighbours are not linearly independent') 

# Determine which pairwise interactions are inferrable
# this is done species by species
Q <- sapply(key_speciesID, function(f){
  
  N_i <- as.matrix(df[df$focal == f, 5:56])
  X_i <- cbind(1,N_i)
  R_i <- pracma::rref(X_i)
  Z_i <- t(R_i) %*% R_i
  
  # param k is inferrable if its corresponding row/column is all 0 except for the k'th element
  # ignore intercept because we always want to include it
  sapply(seq(2, dim(Z_i)[1], 1), function(k){ 
    ifelse(Z_i[k, k] == 1 & sum(Z_i[k, -k]) == 0, 1, 0)
  }) # inferrable params == 1
  
})
Q <- t(Q)
# Q is a matrix of focal x neighbours, if Q[i, j] = 1 then the interaction between i and j is

# transform data into format required by STAN
stan.data <- rem_dataprep(fecundities)
stan.data$Q <- Q

message(paste0('Community selected: ', comm))
message(paste0('Fecundity data dimensions = ', dim(fecundities)[1], ', ', dim(fecundities)[2]))
message(paste0('Number of focals = ', length(key_speciesID)))
message(paste0('Number of neighbours = ', length(key_neighbourID)))

#--------------------------------------------------
# Estimate interactions with a joint NDD*RI model |
#--------------------------------------------------

fit <- stan(file = 'joint_model.stan',
            data =  stan.data,               # named list of data
            chains = 1,
            warmup = 5000,          # number of warmup iterations per chain
            iter = 10000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10)
)

# parameters of interest
param.vec <- c('beta_i0', 'beta_ij', 'effect', 'response', 're', 'inter_mat',
               'mu', 'disp_dev', 'sigma')

# Raw output
#------------
save(fit, file = paste0('model/output/model_fit.Rdata')) # model fit
# Save the 'raw' draws from the posterior
joint.post.draws <- extract.samples(fit)
save(joint.post.draws, file = paste0('model/output/post_draws.Rdata'))

# Save mean, 10% and 90% quantiles for each parameter, as well as n_eff and Rhat
fit_sum <- summary(fit, pars = param.vec, probs = c(0.1, 0.9))$summary
write.csv(fit_sum, file = paste0('model/output/summary_of_draws.csv'), row.names = T)

# Save the logarithm of the (unnormalized) posterior density (lp__)
log_post <- unlist(extract(fit, 'lp__'))
write.csv(log_post, file = paste0('model/output/log_post.csv'), row.names = F)

# Validation
#------------
# Get Geweke statistics
matrix_of_draws <- as.matrix(fit)
gew <- coda::geweke.diag(matrix_of_draws)
write.csv(gew$z, 'model/validation/gew_stats.csv')
# get adjusted values (see boral() package)
gew.pvals <- 2*pnorm(abs(unlist(gew$z)), lower.tail = FALSE)
adj.gew <- p.adjust(gew.pvals, method = "holm")
write.csv(adj.gew, 'model/validation/gew_stats_holmadjust.csv')
print(paste0('Range of p-values for chain convergence: ', min(na.omit(adj.gew)), ' to ',  max(na.omit(adj.gew))))

png('model/validation/geweke_dist.png', width = 500, height = 500)
plot(density(na.omit(gew$z)))
lines(density(rnorm(10000)), col = 'red')
dev.off()
# Diagnostics
stan_diagnostic(fit, 'model/validation/')
# Traceplots and posterior uncertainty intervals
stan_model_check(fit, 'model/validation/', params = param.vec)
# Posterior predictive check
stan_post_pred_check(joint.post.draws, 'model/validation/', stan.data)

# Parameter outputs - draw 1000 samples from the 80% confidence intervals and save 
#------------------
# sample raw model params
sapply(param.vec[param.vec != 'sigma' & param.vec != 'inter_mat'], function(p) {
  
  p.samples <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  })
  write.csv(p.samples, paste0('model/output/', p, '_samples.csv'), 
            row.names = F)
  
})
# sigma alph mucks up because it is one value only, ifm_alphas mucks up because it's a 3d array
sigma <- joint.post.draws$sigma
sigma <- sample(sigma[sigma > quantile(sigma, 0.1) & sigma < quantile(sigma, 0.9)], size = 1000)
write.csv(sigma, paste0('model/output/sigma_samples.csv'), 
          row.names = F)


# Transformed parameters
#-----------------------
# Intrinsic growth rate (lambda)
growth.rates.samples <- apply(joint.post.draws$beta_i0, 2, function(x){
  sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
})
# exponentiate to get lambda
growth.rates.samples <- exp(growth.rates.samples)
write.csv(growth.rates.samples, paste0('model/transformed/lambda_samples.csv'), row.names = F)


# Interactions (alphas)
# get posterior draws for all interactions (from the IFM)
alphas <- joint.post.draws$inter_mat
# columns are interactions (in same order as the interactions vector - focals, then neighbours)
alphas <- as.data.frame(aperm(alphas, perm = c(1, 3, 2)))
colnames(alphas) <- grep('ifm_alpha', rownames(fit_sum), value = T)
# take the 80% posterior interval
alphas <- apply(alphas, 2, function(x) {
  inter <- x[x > quantile(x, 0.1) & x < quantile(x, 0.9)]
  if (length(inter > 0)) {sample(inter, size = 1000)} else {rep(0, 1000)}
  # this is for those unobserved interactions (0)
})

# calculate interactions according to the response effect model
response <- apply(joint.post.draws$response, 2, function(x) {
  sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
})
effect <- apply(joint.post.draws$effect, 2, function(x) {
  sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
})

re_alphas <- lapply(c(1:length(key_speciesID)), function(x) {
  ri <- response[, x]
  ri <- sample(ri, length(ri))  # randomly re-order the ri vector
  return(ri*effect)})
re_alphas <- do.call(cbind, re_alphas)

# Verify!
png(paste0('model/validation/ifm_vs_rem_alphas.png'))
plot(alphas, re_alphas, xlab = 'IFM alphas', ylab='Response*Effect',
     xlim = c(min(re_alphas), max(re_alphas)),
     ylim = c(min(re_alphas), max(re_alphas)))
abline(0,1)
dev.off()

# get unobserved estimates only
unobs <- re_alphas[ , apply(alphas, 2, function(x) {all(x == 0)})]
png(paste0('model/validation/alpha_est_distr.png'))
par(mfrow=c(3,1))
hist(alphas[  , apply(alphas, 2, function(x) {all(x != 0)})], xlab = "", breaks = 30,
     main = "IFM alphas (0's removed)", xlim = c(min(re_alphas), max(re_alphas)))
hist(re_alphas[ , apply(alphas, 2, function(x) {all(x != 0)})],  xlab = "", breaks = 30,
     main = 'REM alphas - Unrealised estimates removed', xlim = c(min(re_alphas), max(re_alphas)))
hist(unobs,  xlab = "", main = 'REM alphas - Unrealised interactions only', breaks = 30,
     xlim = c(min(re_alphas), max(re_alphas)))
dev.off()

# replace unobserved interactions (0 in alphas) with the values predicted by the rem
alphas[ , apply(alphas, 2, function(x) {all(x == 0)})] <- 
  re_alphas[ , apply(alphas, 2, function(x) {all(x == 0)})]
write.csv(alphas, paste0('model/transformed/alpha_samples.csv'), row.names = F)

# Scale the alphas and save
scaled_alphas <- scale_interactions(alphas, growth.rates.samples, key_speciesID, key_neighbourID, comm)
save(scaled_alphas, file = paste0('model/transformed/scaled_alpha_matrices.Rdata')) 

Sys.time()




