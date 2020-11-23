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

# get real data
fecundities <- read.csv('3.case_study/data/fecundities0.csv', stringsAsFactors = F)
seeds <- fecundities$seeds
focalobs <- as.numeric(as.factor(fecundities$focal))

# extract mu and phi
mu <- p.samples$mu
phi <- (p.samples$disp_dev^2)^(-1)

# generating posterior predictions
seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
for (i in 1:dim(mu)[1]) {     # for each posterior draw
  for (j in 1:dim(mu)[2]) {    # for each observation 
    # draw from the predicted distribution
    seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i, focalobs[j]])  
  }
}
# and using just the mean of the parameters
m.mu <- colMeans(mu)
m.phi <- colMeans(phi)
mean_seed_pred <- sapply(1:length(m.mu) , function(x) {
  rnbinom(1, mu = m.mu[x], size = m.phi[focalobs[x]])
})

# get maximum density for plot limits
max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                     max(density(seeds)$y)))


png('2.analyses/figures_mss/postpredch_nolimit.png', width=800, height=800)
# start a plot with the first draw 
ppc.plot <- plot(density(seed_pred[1, ]), 
               #  xlim = c(0, 1000), # this is only so we can zoom in
                 ylim = c(0, max.density), 
                 col = 'lightgrey',
                 ylab = 'Seed density',
                 main = 'Post. pred. check',
                 sub = '(grey = predicted, black = observed)') 
for (i in 2:dim(seed_pred)[1]) {
  # add a line for each draw
  ppc.plot <- lines(density(seed_pred[i, ]), col = 'lightgrey')
}

polygon(c(density(seed_pred)$x, density(x2)$x),  # X-Coordinates of polygon
        c(density(seed_pred)$y, density(x2)$y),  # Y-Coordinates of polygon
        col = "#1b98e0") 


# add the 'mean prediction'
ppc.plot <- lines(density(mean_seed_pred), col = 'grey40', lwd = 1.5)  
# add the actual data
ppc.plot <- lines(density(seeds), col = 'black')  
print(ppc.plot)
dev.off()



#-------------------------------
# Figure 3: Response vs Impact |
#-------------------------------
S <- dim(p.samples$a)[[2]]

png('2.analyses/figures_mss/response_impact2.png', width = 400, height = 700)
plot(p.samples$response, p.samples$effect[ , 1:S], type = 'n',
     xlab = expression(italic('response i')), ylab = expression(italic('impact i')))
abline(h = 0, lty = 2)
points(p.samples$response, p.samples$effect[ , 1:S],
     pch = 16, 
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2))
points(colMeans(p.samples$response), colMeans(p.samples$effect[ , 1:S]),
       pch = 18, cex = 2)
dev.off()

#--------------------------------------
# Figure 4: Pretty network vs cooccur |
#--------------------------------------





#-------------------------------------------
# Figure 5: Applications - species effects |
#-------------------------------------------



