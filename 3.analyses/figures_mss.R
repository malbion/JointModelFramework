# CompNet Methods paper
# Script to create figures for manuscript 

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')


# Figures 1-3 use 'raw' parameters extracted from the model output 
# Figures 4 and 5 used interaction estimates which have been rescaled with a pop dyn model


############################################################################################
#                               PREP FOR FIGURES 1-3
############################################################################################

# # Get raw quantities from model output
# load('2.case_study/model/output/post_draws.Rdata')
# 
# # 80% posterior of relevant parameters - vectors
# param.vec <- c('beta_i0', 'beta_ij', 'effect', 'response', 're', 'mu', 'disp_dev')
# p.samples <- list()
# p.samples <- sapply(param.vec, function(p) {
#   p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
#     sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
#   }) 
# })

rim_mat <- as.matrix(read.csv('2.case_study/model/output/RIM_betaij_samples.csv'))
nddm_mat <- as.matrix(read.csv('2.case_study/model/output/NDDM_betaij_samples.csv'))


#---------------------------------
# Figure 1: IFM vs RIM estimates |
#---------------------------------

png('3.analyses/figures_mss/interaction_estimates.png', width = 1400, height = 700)
par(mfrow=c(2,2), cex = 1.5)
hist(nddm_mat[ , apply(nddm_mat, 2, function(x) {all(x != 0)})], 
     xlab = "", breaks = 30,
     main = "Realised interactions (NDDM)", 
     xlim = c(min(rim_mat), max(rim_mat)))
plot(rim_mat[nddm_mat!=0], nddm_mat[nddm_mat!=0], 
     ylab = 'NDDM estimates', xlab='Response*Impact estimates',
     xlim = c(min(rim_mat), max(rim_mat)), ylim = c(min(rim_mat), max(rim_mat)), 
     pch = 16, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
abline(0,1)
hist(rim_mat[ , apply(nddm_mat, 2, function(x) {all(x == 0)})],  
     xlab = "", breaks = 30,
     main = 'Unrealised interactions (RIM)', 
     xlim = c(min(rim_mat), max(rim_mat)))
hist(rim_mat[ , apply(nddm_mat, 2, function(x) {all(x != 0)})],  
     xlab = "", breaks = 30,
     main = 'Realised interactions (RIM)', 
     xlim = c(min(rim_mat), max(rim_mat)))
dev.off()


#---------------------------------------
# Figure 2: Posterior predictive check |
#---------------------------------------

# get real data
fecundities <- read.csv('2.case_study/data/fecundities0.csv', stringsAsFactors = F)
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
# log transform seed predictions
seed_pred <- log(seed_pred)
# and using just the mean of the parameters
# 06/2021 changed to medians
# m.mu <- colMeans(mu)
# m.phi <- colMeans(phi)
m.mu <- apply(mu, 2, median)
m.phi <- apply(phi, 2, median)
mean_seed_pred <- sapply(1:length(m.mu) , function(x) {
  rnbinom(1, mu = m.mu[x], size = m.phi[focalobs[x]])
})
mean_seed_pred <- log(mean_seed_pred)
# # get the seed density for each sample 
# seed_dens <- apply(seed_pred, 1, function(x) {list(density(x)$x, density(x)$y)})
# 
# seed_densX <- apply(seed_pred, 1, function(x) {density(x)$x})
# seed_densY <- apply(seed_pred, 1, function(x) {density(x)$y})
# # get the minimum and maximum y for each row
# minY <- apply(seed_densY, 1, quantile, 0.025)
# maxY <- apply(seed_densY, 1, quantile, 0.975)

# sapply(1:nrow(seed_densY), function(x) {
#   temp1 <- quantile(seed_densY[x, ], 0.025)
#   temp2 <- seed_densX[x, ]
# })


# get maximum density for plot limits
max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                     max(density(seeds)$y)))


png('3.analyses/figures_mss/postpredch.png', width=500, height=500)
# start a plot with the first draw 
ppc.plot <- plot(density(seed_pred[1, ]), 
                 type = 'n',
                 bty = 'n',
                 ylim = c(0, max.density+0.03), 
                 col = 'lightgrey',
                 ylab = 'Seed probability density',
                 xlab = 'Log seed production',
                 main = '',
                 # sub = '(grey = predicted, red = observed)'
                 ) 
for (i in 1:dim(seed_pred)[1]) {
  # add a line for each draw
  ppc.plot <- lines(density(seed_pred[i, ]), col = 'lightgrey')
}
# add the 'mean prediction'
ppc.plot <- lines(density(mean_seed_pred), col = 'black', lwd = 1)  
# add the actual data
ppc.plot <- lines(density(log(seeds)), col = 'red', lwd = 1)  
# add legend
ppc.plot <- legend('topright', 
                  legend = c('predicted', 'observed'), 
                  col = c('lightgrey', 'red'), lwd = c(6, 1),
                  bty = 'n')
print(ppc.plot)
dev.off()


#-------------------------------
# Figure 3: Response vs Impact |
#-------------------------------

S <- dim(p.samples$beta_i0)[[2]]

png('3.analyses/figures_mss/response_impact2.png', width = 400, height = 700)
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


load('2.case_study/model/transformed/scaled_betaij_matrices.Rdata')

#--------------------------------------
# Figure 4: Pretty network vs cooccur |
#--------------------------------------

# Median interaction estimates for joint model
bij.med <- apply(scaled_betas, c(1, 2), median)
bij.med <- bij.med[ , 1:nrow(bij.med)]

# set up colours for nodes (linking to next figure)
invasives <- c('ARCA', 'PEAI', 'HYPO')
foundation <- c('VERO', 'POCA')
keyst <- c('GITE', 'TROR', 'HAOD')
all.sp <- rep('white', length(dimnames(bij.med)$species))
names(all.sp) <- dimnames(bij.med)$species
# all.sp[invasives] <- 'tomato'
all.sp[foundation] <- 'plum'
all.sp[keyst] <- 'seagreen2'
# node.outline.width <- rep(1, length(dimnames(bij.med)$species))
# names(node.outline.width) <- dimnames(bij.med)$species
# node.outline.width[ c(invasives, foundation, keyst)] <- 3

# get the cooccurence matrix
cooc <- read.csv('3.analyses/0_cooc.csv', stringsAsFactors = F, row.names = 1)
cooc <- cooc[1:nrow(bij.med), 1:nrow(bij.med)]

library(qgraph)

png('3.analyses/figures_mss/networks_strengths_med.png', 
    width = 600, height = 800, units = 'px')
par(mfrow=c(2, 1))
# plot all interactions
qgraph(bij.med,  # plot interaction means
     #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
       layout = 'circle',
       negCol = 'royalblue4',   # facilitation = blue
       posCol = 'orange',       # competition = orange
       color = all.sp,
       labels = dimnames(bij.med)$species,
       fade = T, directed = T,
       title = 'A', title.cex =5)
# plot from the cooccur package
qgraph(cooc,
       layout = 'circle', 
       negCol = 'orange',   # swap the colours around
       posCol = 'royalblue4',     
       color = all.sp,
       labels = dimnames(bij.med)$species,
       fade = T,
       title = 'B', title.cex =5)
dev.off()

png('3.analyses/figures_mss/networks_C_F_cooc.png', 
    width = 600, height = 1200, units = 'px')
par(mfrow=c(3, 1))
# plot competition only
qgraph(bij.med,  # plot interaction means
       #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
       layout = 'circle',
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),   # facilitation = transaprent
       posCol = 'orange',       # competition = orange
       color = all.sp,
      # border.width = node.outline.width,
       labels = dimnames(bij.med)$species,
       fade = T, directed = T,
       title = 'A', title.cex =5)
# plot facilitation only
qgraph(bij.med,  # plot interaction means
       #  edge.width = (alpha_var*100),  # set edge width to be equal to the variance
       layout = 'circle',
       negCol = 'royalblue4',   # facilitation = transaprent
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),       # competition = orange
       color = all.sp,
       labels = dimnames(bij.med)$species,
       fade = T, directed = T,
       title = 'B', title.cex =5)
# plot from the cooccur package
qgraph(cooc,
       layout = 'circle', 
       negCol = 'orange',   # swap the colours around
       posCol = 'royalblue4',     
       color = all.sp,
       labels = dimnames(bij.med)$species,
       fade = T,
       title = 'C', title.cex =5)
dev.off()

#-------------------------------------------
# Figure 5: Applications - species effects |
#-------------------------------------------

sp.abunds <- read.csv('2.case_study/data/plot_species_abundances.csv', stringsAsFactors = F)
sp.abunds <- split(sp.abunds, as.factor(sp.abunds$species))
sp.abunds <- do.call(rbind,
                     lapply(sp.abunds, function(x) sum(x[ ,'Num_indivs'])))
sp.abunds <- sp.abunds[rownames(bij.med), ]
# log scale 
sp.abunds <- log(sp.abunds)
sp.abunds.rep <- matrix(data = rep(sp.abunds, dim(scaled_betas)[3]), 
                        nrow = length(sp.abunds), ncol = dim(scaled_betas)[3])

# get the effects of species on themselves 
intras <- apply(scaled_betas, 3, diag)
upper.intra <- apply(intras, 1, quantile, 0.75)
lower.intra <- apply(intras, 1, quantile, 0.25)
mean.intra <- apply(intras, 1, median)

# get the effects of species on others
inters.out <- scaled_betas[ , 1:nrow(scaled_betas), ]
sum.out <- apply(inters.out, c(2, 3), sum)
sum.out <- sum.out - intras   # remove intraspecific interactions
upper.sum.out <- apply(sum.out, 1, quantile, 0.75)
lower.sum.out <- apply(sum.out, 1, quantile, 0.25)
mean.sum.out <- apply(sum.out, 1, median)

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
# median values for each species
med.comp <- apply(sum.comp, 1, median)
med.faci <- apply(sum.faci, 1, median)

# get 'special' species
invasives <- c('ARCA', 'PEAI', 'HYPO')
foundation <- c('VERO', 'POCA')
keyst <- c('GITE',  'HAOD')

###### PLOT
png('3.analyses/figures_mss/species_effects.png', 
    width = 500, height = 1500, units = 'px')
par(mfrow=c(3,1), cex=1.2)

# 1. intra vs abundance
#-------------------
plot(intras, sp.abunds.rep,
     xlab = 'Intraspecific interactions (self-regulation)',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2)
title(main = 'A', adj = 0)
polygon(x = c(-0.129 , 0, 0, -0.129), y = c(4.4, 4.4, 9.3, 9.3),
        col = 'ivory2', border = NA)
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
       bg = 'black', cex = 1)
points(mean.intra[foundation], sp.abunds[foundation], pch = 24,
       bg = 'purple', cex = 1.8)
# points(mean.intra[keyst], sp.abunds[keyst], pch = 23, 
#        bg = 'orange', cex = 1.3)
# points(mean.intra[invasives], sp.abunds[invasives], pch = 23, 
#        bg = 'red', cex = 1.3)
text(mean.intra[foundation], sp.abunds[foundation], 
     labels = names(mean.intra[foundation]), 
     pos = 4, col = 'darkorchid4', offset = 1, cex = 1.4)


# 2. out-strength vs abundance 
#-----------------------------
plot(sum.out, sp.abunds.rep, 
     xlab = 'Net effect on neighbours (Out-strength)',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2)
title(main = 'B', adj = 0)
polygon(x = c(-0.475 , 0, 0, -0.475), y = c(4.4, 4.4, 9.3, 9.3),
        col = 'ivory2', border = NA)
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
       bg = 'black', cex = 1)
# points(mean.sum.out[foundation], sp.abunds[foundation], pch = 23, 
#        bg = 'royalblue', cex = 1.3)
points(mean.sum.out[keyst], sp.abunds[keyst], pch = 24, 
       bg = 'chartreuse3', cex = 1.8)
# ppoints(mean.sum.out[invasives], sp.abunds[invasives], pch = 23, 
#        bg = 'red', cex = 1.3)
lab.x.pos <- mean.sum.out[keyst]
# lab.x.pos['TROR'] <- lab.x.pos['TROR'] - 0.15
text(lab.x.pos, sp.abunds[keyst], 
     labels = names(mean.sum.out[keyst]), pos = 3, col = 'green4', 
     cex = 1.4, offset = 0.7)


# 3. competitive vs facilitative effects
#-------------------------------------
plot(sum.comp, -sum.faci,  
     las = 1, bty = 'n', cex.lab = 1.2,
     # type = 'n',
     pch = 16, cex=0.7, col = 'grey',
     xlab = 'Absolute sum of competitive effects', 
     ylab = 'Absolute sum of facilitative effects')
title(main = 'C', adj = 0)
abline(v=median(sum.comp), lty = 2)
abline(h=median(-sum.faci), lty = 2)
# points(sum.comp, -sum.faci, las = 1, bty = 'n', 
#        pch = 16, cex=0.7, col = 'grey')
points(med.comp, -med.faci, 
       pch = 23, cex=1, bg='black')
# points(rowMeans(sum.comp)[foundation], -rowMeans(sum.faci)[foundation], 
#        pch = 23, bg = 'royalblue', cex = 1.3)
# points(rowMeans(sum.comp)[keyst], -rowMeans(sum.faci)[keyst], 
#        pch = 23, bg = 'orange', cex = 1.3)
points(med.comp[invasives], -med.faci[invasives],  
       pch = 24, bg = 'red', cex = 1.8)
lab.x.pos <- med.comp[invasives] 
lab.x.pos['HYPO'] <- lab.x.pos['HYPO'] + 0.1
lab.x.pos['PEIA'] <- lab.x.pos['PEIA'] - 0.3
# lab.x.pos['ARCA'] <- lab.x.pos['ARCA'] + 0.1
lab.y.pos <- (-med.faci[invasives] + 0.065)
lab.y.pos['PEAI'] <- lab.y.pos['PEAI'] +0.05
text(lab.x.pos, lab.y.pos, 
     labels = names(med.comp[invasives]), 
     pos =  4, col = 'darkred', offset = 0.3, cex = 1.4)

dev.off()


# Another figure - cooccur strengths vs log abundance? 
#------------------------------------------------------
png('3.analyses/figures_mss/cooccur_vs_abund.png', 
    width = 500, height = 500, units = 'px')
par(cex=1.2)
plot(-colSums(cooc, na.rm = T), sp.abunds, 
     xlab = 'Absolute sum of association strengths',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2,
     xlim = c(0, 40))
abline(v=median(-colSums(cooc, na.rm = T)), lty = 2)
abline(h=median(sp.abunds), lty = 2)
# points for species means
points(-colSums(cooc, na.rm = T), sp.abunds, pch = 23, 
       bg = 'black', cex = 1)
# points(mean.sum.out[foundation], sp.abunds[foundation], pch = 23, 
#        bg = 'royalblue', cex = 1.3)
points(-colSums(cooc, na.rm = T)[keyst], sp.abunds[keyst], pch = 24, 
       bg = 'chartreuse3', cex = 1.8)
# ppoints(mean.sum.out[invasives], sp.abunds[invasives], pch = 23, 
#        bg = 'red', cex = 1.3)
text(-colSums(cooc, na.rm = T)[keyst], sp.abunds[keyst], 
     labels = names(-colSums(cooc, na.rm = T)[keyst]), pos = 4, col = 'chartreuse4', 
     cex = 1.4, offset = 1.7)
dev.off()

