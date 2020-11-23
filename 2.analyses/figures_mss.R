# CompNet Methods paper
# Script to create figures for manuscript 

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')


# Figures 1-3 use 'raw' parameters extracted from the model output 
# Figures 4 and 5 used interaction estimates which have been rescaled with a pop dyn model

# 1. Get raw quantities from model output
#-----------------------------------------
load('../2018_Compnet/stormland/model/output/0/post_draws.Rdata')

# 80% posterior of relevant parameters - vectors
param.vec <- c('a', 'interactions', 'effect', 'response', 're', 'mu', 'disp_dev')
p.samples <- list()
p.samples <- sapply(param.vec, function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  }) 
})
# get ALL interaction estimates for IFM 
ifm_mat <- as.data.frame(aperm(ifm_mat, perm = c(1, 3, 2)))
# take the 80% posterior interval
ifm_mat <- apply(ifm_mat, 2, function(x) {
  inter <- x[x > quantile(x, 0.1) & x < quantile(x, 0.9)]
  # this is for those unobserved interactions (0):
  if (length(inter > 0)) {sample(inter, size = 1000)} else {rep(0, 1000)} 
})
# get ALL interaction estimates for REM 
rim_mat <- lapply(c(1:dim(p.samples[['a']])[[2]]), function(x) {
  resp <- p.samples$response[, x]
  resp <- sample(resp, length(resp))  # randomly re-order the response vector
  return(resp*p.samples$effect)})
rim_mat <- do.call(cbind, rim_mat)




#---------------------------------
# Figure 1: IFM vs RIM estimates |
#---------------------------------

png('2.analyses/figures_mss/interaction_estimates.png', width = 1400, height = 700)
par(mfrow=c(2,2), cex = 1.5)
hist(ifm_mat[ , apply(ifm_mat, 2, function(x) {all(x != 0)})], 
     xlab = "", breaks = 30,
     main = "Observed interactions (IFM)", 
     xlim = c(min(rim_mat), max(rim_mat)))
plot(rim_mat[ifm_mat!=0], ifm_mat[ifm_mat!=0], 
     ylab = 'IFM betas', xlab='Response*Effect estimates',
     xlim = c(min(rim_mat), max(rim_mat)), ylim = c(min(rim_mat), max(rim_mat)), 
     pch = 16, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
abline(0,1)
hist(rim_mat[ , apply(ifm_mat, 2, function(x) {all(x == 0)})],  
     xlab = "", breaks = 30,
     main = 'Unrealised interactions (REM)', 
     xlim = c(min(rim_mat), max(rim_mat)))
hist(rim_mat[ , apply(ifm_mat, 2, function(x) {all(x != 0)})],  
     xlab = "", breaks = 30,
     main = 'Observed interactions (REM)', 
     xlim = c(min(rim_mat), max(rim_mat)))
dev.off()



#---------------------------------------
# Figure 2: Posterior predictive check |
#---------------------------------------




#-------------------------------
# Figure 3: Response vs Impact |
#-------------------------------





#--------------------------------------
# Figure 4: Pretty network vs cooccur |
#--------------------------------------





#-------------------------------------------
# Figure 5: Applications - species effects |
#-------------------------------------------



