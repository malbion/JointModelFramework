# Objective: running a joint model combining NDDM and RIM to estimate all interactions at the same time. 
# 
# Run in Terminal: 
#   
#   R CMD BATCH '--args 0' model.R model/0_model_script.Rout
# 
# substituting '0' for the comm number (0, 1, 2, or 3) 
# 
# This calls up model.R, which runs joint_model.stan. 
# results got to model/output, /transformed and /validation


# Run the joint model: 
# 1. Estimate inferrable interactions with NDDM 
# 2. Estimate response and effect interactions when non-inferrable

# Then: 
# 3. Extract the intrinsic growth rates 
# 5. Scale the interactions into per-capita growth rates

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

setwd(here('2.case_study/'))

source('../1.code/data_prep.R')
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

# transform data into format required by STAN
stan.data <- data_prep(perform = 'seeds', focal = 'focal', 
                       nonNcols = 4, df = fecundities)


message(paste0('Community selected: ', comm))
message(paste0('Fecundity data dimensions = ', dim(fecundities)[1], ', ', dim(fecundities)[2]))
message(paste0('Number of focals = ', length(key_speciesID)))
message(paste0('Number of neighbours = ', length(key_neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$K)))

#--------------------------------------------------
# Estimate interactions with a joint NDD*RI model |
#--------------------------------------------------

fit <- stan(file = '../1.code/joint_model.stan',
            data =  stan.data,               # named list of data
            chains = 1,
            warmup = 5000,          # number of warmup iterations per chain
            iter = 10000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10)
)

# parameters of interest
param.vec <- fit@model_pars[!fit@model_pars %in% c('response1', 'responseSm1', 'lp__')]


# Raw output
#------------
save(fit, file = paste0('model/output/model_fit.Rdata')) # model fit
# Save the 'raw' draws from the posterior
post.draws <- extract.samples(fit)
save(post.draws, file = paste0('model/output/post_draws.Rdata'))

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
abline(v = -2, lty = 2)
abline(v = 2, lty = 2)
dev.off()


# Diagnostics
stan_diagnostic(fit, 'model/validation/')
# Traceplots and posterior uncertainty intervals
stan_model_check(fit, 'model/validation/', params = param.vec)
# Posterior predictive check
stan_post_pred_check(post.draws, 'mu', 'model/validation/', stan.data)
stan_post_pred_check(post.draws, 'mu2', 'model/validation/', stan.data)

# Parameter outputs - draw 1000 samples from the 80% confidence intervals and save 
#------------------
# this works for parameters that are not the interaction matrices
sapply(param.vec[!param.vec %in% c('ri_betaij', 'ndd_betaij')], function(p) {
  
  p.samples <- apply(post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  })
  write.csv(p.samples, paste0('model/output/', p, '_samples.csv'), 
            row.names = F)
  
})

# Interactions (betaij)
# joint interactions 
betaij <- apply(post.draws$ndd_betaij, c(2, 3), function(x) {
  sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
})
betaijS <- as.data.frame(aperm(betaij, perm = c(1, 3, 2)))
colnames(betaijS) <- grep('beta_ij', rownames(fit_sum), value = T)
write.csv(betaijS, paste0('model/output/joint_betaij_samples.csv'), row.names = F)

# rim interactions only
rim_betaij <- apply(post.draws$ri_betaij, c(2, 3), function(x) {
  sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
})
rim_betaijS <- as.data.frame(aperm(rim_betaij, perm = c(1, 3, 2)))
colnames(rim_betaijS) <- grep('ri_betaij', rownames(fit_sum), value = T)
write.csv(rim_betaijS, paste0('model/output/RIM_betaij_samples.csv'), row.names = F)

# replace uninferrable interactions with 0 to get the NDDM estimates only 
nddm_betaij <- as.data.frame(t(apply(betaij, 1, function(x) x*stan.data$Q)))  
colnames(nddm_betaij) <- grep('ndd_betaij', rownames(fit_sum), value = T)
write.csv(nddm_betaij, paste0('model/output/NDDM_betaij_samples.csv'), row.names = F)

# Let's do a few plots while we're at it
# interactions inferred by the NDDM only
nddm_betaij_inf <-as.matrix(nddm_betaij[  , which(colSums(nddm_betaij) != 0)]) # NDDM without non-inferrables
# RIM estimates for those inferrable interactions
rim_betaij_inf <- as.matrix(rim_betaijS[  , which(colSums(nddm_betaij) != 0)]) # NDDM without non-inferrables
# RIM estimates for non-inferrable interactions only
rim_betaij_noinf <- as.matrix(rim_betaijS[  , which(colSums(nddm_betaij) == 0)])

# Check estimates of inferrable interactions from both models 
png(paste0('model/validation/nddm_vs_rim_alphas.png'))
plot(nddm_betaij_inf, rim_betaij_inf, 
     xlab = 'NDDM interactions (inferrable only)', 
     ylab = 'RIM interactions (inferrable only)',
     xlim = c(min(nddm_betaij), max(nddm_betaij)),
     ylim = c(min(nddm_betaij), max(nddm_betaij)))
abline(0,1)
dev.off()

# check distribution of inferrable and non-inferrable interactions
png(paste0('model/validation/betaij_est_distr.png'))
par(mfrow=c(3,1))
hist(nddm_betaij_inf, xlab = "", breaks = 30,
     main = "Inferrable interactions (NDDM)", xlim = c(min(betaijS), max(betaijS)))
hist(rim_betaij_inf,  xlab = "", breaks = 30,
     main = 'Inferrable interactions (RIM)', xlim = c(min(betaijS), max(betaijS)))
hist(rim_betaij_noinf,  xlab = "", breaks = 30,
     main = 'Non-inferrable interactions (RIM)', xlim = c(min(betaijS), max(betaijS)))
dev.off()

# Transformed parameters
#-----------------------
# Intrinsic growth rate (lambda)
growth.rates.samples <- apply(post.draws$beta_i0, 2, function(x){
  sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
})
# exponentiate to get lambda
growth.rates.samples <- exp(growth.rates.samples)
write.csv(growth.rates.samples, paste0('model/transformed/lambda_samples.csv'), row.names = F)

# Scale the alphas and save
scaled_betas <- scale_interactions(betaijS, growth.rates.samples, key_speciesID, key_neighbourID, comm)
save(scaled_betas, file = paste0('model/transformed/scaled_betaij_matrices.Rdata')) 

Sys.time()




