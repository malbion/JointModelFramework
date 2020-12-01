# CompNet Methods paper
# Script to create figures for manuscript 

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')


# Figures 1-3 use 'raw' parameters extracted from the model output 
# Figures 4 and 5 used interaction estimates which have been rescaled with a pop dyn model


############################################################################################
#                               PREP FOR FIGURES 1-3
############################################################################################

# Get raw quantities from model output
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
ifm_mat <- as.data.frame(aperm(joint.post.draws$ifm_alpha, perm = c(1, 3, 2)))
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
     main = "Observed interactions (NDDM)", 
     xlim = c(min(rim_mat), max(rim_mat)))
plot(rim_mat[ifm_mat!=0], ifm_mat[ifm_mat!=0], 
     ylab = 'NDDM betas', xlab='Response*Impact estimates',
     xlim = c(min(rim_mat), max(rim_mat)), ylim = c(min(rim_mat), max(rim_mat)), 
     pch = 16, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
abline(0,1)
hist(rim_mat[ , apply(ifm_mat, 2, function(x) {all(x == 0)})],  
     xlab = "", breaks = 30,
     main = 'Unrealised interactions (RIM)', 
     xlim = c(min(rim_mat), max(rim_mat)))
hist(rim_mat[ , apply(ifm_mat, 2, function(x) {all(x != 0)})],  
     xlab = "", breaks = 30,
     main = 'Observed interactions (RIM)', 
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

# get the seed density for each sample 
seed_dens <- apply(seed_pred, 1, function(x) {list(density(x)$x, density(x)$y)})

seed_densX <- apply(seed_pred, 1, function(x) {density(x)$x})
seed_densY <- apply(seed_pred, 1, function(x) {density(x)$y})
# get the minimum and maximum y for each row
minY <- apply(seed_densY, 1, quantile, 0.025)
maxY <- apply(seed_densY, 1, quantile, 0.975)

sapply(1:nrow(seed_densY), function(x) {
  temp1 <- quantile(seed_densY[x, ], 0.025)
  temp2 <- seed_densX[x, ]
})


# get maximum density for plot limits
max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                     max(density(seeds)$y)))


png('2.analyses/figures_mss/postpredch.png', width=500, height=500)
# start a plot with the first draw 
ppc.plot <- plot(density(seed_pred[1, ]), 
                 type = 'n',
                  xlim = c(0, 400), # this is only so we can zoom in
                 ylim = c(0, max.density), 
                 col = 'lightgrey',
                 ylab = 'Seed density',
                 main = 'Post. pred. check',
                 sub = '(grey = predicted, black = observed)') 
for (i in 1:dim(seed_pred)[1]) {
  # add a line for each draw
  ppc.plot <- lines(density(seed_pred[i, ]), col = 'lightgrey')
}
# ppc.plot <- points(seed_densX, seed_densY, pch = 16, 
#                    col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
# polygon(c(meanX, rev(meanX)),  # X-Coordinates of polygon
#         c(minY, rev(maxY)),  # Y-Coordinates of polygon
#         col = "grey") 


# add the 'mean prediction'
ppc.plot <- lines(density(mean_seed_pred), col = 'black', lwd = 1)  
# add the actual data
ppc.plot <- lines(density(seeds), col = 'red', lwd = 1)  
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




rm(joint.post.draws)


############################################################################################
#                               PREP FOR FIGURES 4 & 5
############################################################################################


load('3.case_study/model/transformed/scaled_alpha_matrices.Rdata')



#--------------------------------------
# Figure 4: Pretty network vs cooccur |
#--------------------------------------

# get the mean and variance estimates for interactions
alpha_means <- apply(scaled_alphas, c(1, 2), mean)
alpha_means <- alpha_means[ , 1:nrow(alpha_means)]
alpha_var <- apply(scaled_alphas, c(1, 2), var)
alpha_var <- alpha_var[ , 1:nrow(alpha_var)]

# getting the RIM alphas and variances only
ifm_means <- colMeans(ifm_mat)
ifm_means <- matrix(data = ifm_means, 
                    nrow = dim(p.samples$response)[[2]],
                    ncol = dim(p.samples$effect)[[2]],
                    byrow = T)
ifm_means <- ifm_means[ , 1:nrow(ifm_means)]

alpha_rim <- alpha_means
alpha_rim[which(ifm_means != 0, arr.ind = T)] <- 0
alpha_var_rim <- alpha_var
alpha_var_rim[which(ifm_means != 0, arr.ind = T)] <- 0


# get the cooccurence matrix
cooc <- read.csv('2.analyses/0_cooc.csv', stringsAsFactors = F, row.names = 1)
cooc <- cooc[1:nrow(alpha_means), 1:nrow(alpha_means)]

library(qgraph)

png('2.analyses/figures_mss/networks_strengths.png', 
    width = 600, height = 800, units = 'px')
par(mfrow=c(2, 1))
# plot all interactions
qgraph(alpha_means,  # plot interaction means
     #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
       layout = 'circle',
       negCol = 'royalblue4',   # facilitation = blue
       posCol = 'orange',       # competition = orange
       fade = T, directed = T,
       title = 'A', title.cex =5)
# # plot those from the RIM only 
# qgraph(alpha_rim,  # plot interaction means
#        edge.width = (alpha_var_rim*100),  # set edge width to be equal to the variance
#        layout = 'circle',
#        negCol = 'royalblue4',   # facilitation = blue
#        posCol = 'orange',       # competition = orange
#        fade = T, directed = T,
#        title = 'B', title.cex =5)
# plot from the cooccur package
qgraph(cooc,
       layout = 'circle', 
       negCol = 'orange',   # swap the colours around
       posCol = 'royalblue4',     
       fade = T,
       title = 'B', title.cex =5)
dev.off()

png('2.analyses/figures_mss/networks_C_F_cooc.png', 
    width = 600, height = 1200, units = 'px')
par(mfrow=c(3, 1))
# plot competition only
qgraph(alpha_means,  # plot interaction means
       #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
       layout = 'circle',
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),   # facilitation = transaprent
       posCol = 'orange',       # competition = orange
       fade = T, directed = T,
       title = 'A', title.cex =5)
# plot facilitation only
qgraph(alpha_means,  # plot interaction means
       #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
       layout = 'circle',
       negCol = 'royalblue4',   # facilitation = transaprent
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),       # competition = orange
       fade = T, directed = T,
       title = 'B', title.cex =5)
# plot from the cooccur package
qgraph(cooc,
       layout = 'circle', 
       negCol = 'orange',   # swap the colours around
       posCol = 'royalblue4',     
       fade = T,
       title = 'B', title.cex =5)
dev.off()

#-------------------------------------------
# Figure 5: Applications - species effects |
#-------------------------------------------

sp.abunds <- read.csv('3.case_study/data/plot_species_abundances.csv', stringsAsFactors = F)
sp.abunds <- split(sp.abunds, as.factor(sp.abunds$species))
sp.abunds <- do.call(rbind,
                     lapply(sp.abunds, function(x) sum(x[ ,'Num_indivs'])))
sp.abunds <- sp.abunds[rownames(alpha_means), ]
# log scale 
sp.abunds <- log(sp.abunds)
sp.abunds.rep <- matrix(data = rep(sp.abunds, dim(scaled_alphas)[3]), 
                        nrow = length(sp.abunds), ncol = dim(scaled_alphas)[3])

# get the effects of species on themselves 
intras <- apply(scaled_alphas, 3, diag)
upper.intra <- apply(intras, 1, quantile, 0.75)
lower.intra <- apply(intras, 1, quantile, 0.25)
mean.intra <- apply(intras, 1, mean)

# get the effects of species on others
inters.out <- scaled_alphas[ , 1:nrow(scaled_alphas), ]
sum.out <- apply(inters.out, c(2, 3), sum)
sum.out <- sum.out - intras   # remove intraspecific interactions
upper.sum.out <- apply(sum.out, 1, quantile, 0.75)
lower.sum.out <- apply(sum.out, 1, quantile, 0.25)
mean.sum.out <- apply(sum.out, 1, mean)

# competitive effects only
sum.comp <- apply(inters.out, c(2, 3), function(x) sum(x[x>0]))
# remove competitive intraspecific interaction
compintra <- intras
compintra[compintra<0] <- 0
sum.comp <- sum.comp - compintra   # remove intraspecific interactions
# facilitative effects only 
sum.faci <- apply(inters.out, c(2, 3), function(x) sum(x[x<0]))
# remove facilitative intraspecies interactions
faciintra <- intras
faciintra[faciintra>0] <- 0
sum.faci <- sum.faci - faciintra

# get 'special' species
invasives <- c('ARCA', 'PEAI', 'HYPO')
foundation <- c('VERO', 'POCA')
keyst <- c('GITE', 'TROR', 'HAOD')

# PLOT
png('2.analyses/figures_mss/species_effects.png', 
    width = 1200, height = 400, units = 'px')
par(mfrow=c(1, 3), cex=1.2)

# 1. intra vs abundance
#-------------------
plot(intras, sp.abunds.rep,
     xlab = 'Intraspecific interactions (self-regulation)',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2)
abline(v=median(intras), lty = 2)
abline(h=median(sp.abunds), lty = 2)
points(intras, sp.abunds.rep,
       pch = 16, col = 'grey', cex = 1.5)
# lines for the 50% CI
sapply(1:length(sp.abunds), function(x) {
  lines(x = c(lower.intra[x], upper.intra[x]), 
        y = c(sp.abunds[x], sp.abunds[x]),
        lwd = 1.5)
}) 
# points for species means
points(mean.intra, sp.abunds, pch = 23, 
       bg = 'black', cex = 1.3)
points(mean.intra[foundation], sp.abunds[foundation], pch = 23, 
       bg = 'royalblue', cex = 1.3)
points(mean.intra[keyst], sp.abunds[keyst], pch = 23, 
       bg = 'orange', cex = 1.3)
points(mean.intra[invasives], sp.abunds[invasives], pch = 23, 
       bg = 'red', cex = 1.3)



# 2. out-strength vs abundance 
#-----------------------------
plot(sum.out, sp.abunds.rep,
     xlab = 'Net effect on neighbours (Out-strength)',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2)
abline(v=median(sum.out), lty = 2)
abline(h=median(sp.abunds), lty = 2)
points(sum.out, sp.abunds.rep,
       pch = 16, col = 'grey', cex = 1.5)
# lines for the 50% CI
sapply(1:length(sp.abunds), function(x) {
  lines(x = c(lower.sum.out[x], upper.sum.out[x]), 
        y = c(sp.abunds[x], sp.abunds[x]),
        lwd = 1.5)
}) 
# points for species means
points(mean.sum.out, sp.abunds, pch = 23, 
       bg = 'black', cex = 1.3)
points(mean.sum.out[foundation], sp.abunds[foundation], pch = 23, 
       bg = 'royalblue', cex = 1.3)
points(mean.sum.out[keyst], sp.abunds[keyst], pch = 23, 
       bg = 'orange', cex = 1.3)
points(mean.sum.out[invasives], sp.abunds[invasives], pch = 23, 
       bg = 'red', cex = 1.3)


# 3. competitive vs facilitative effects
#-------------------------------------
plot(sum.comp, -sum.faci, las = 1, bty = 'n', 
     pch = 16, cex=0.7, col = 'grey',
     xlab = 'Sum of competitive effects', 
     ylab = 'Sum of facilitative effects')
points(rowMeans(sum.comp), -rowMeans(sum.faci), 
       pch = 23, cex=1.1, bg='black')
points(rowMeans(sum.comp)[foundation], -rowMeans(sum.faci)[foundation], 
       pch = 23, bg = 'royalblue', cex = 1.3)
points(rowMeans(sum.comp)[keyst], -rowMeans(sum.faci)[keyst], 
       pch = 23, bg = 'orange', cex = 1.3)
points(rowMeans(sum.comp)[invasives], -rowMeans(sum.faci)[invasives],  
       pch = 23, bg = 'red', cex = 1.3)


dev.off()

