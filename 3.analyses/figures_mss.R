# Methods paper
# Script to create figures for manuscript 

library(here)
setwd(here())

# Figures 1 and Supp Figs 1-4 use 'raw' parameters extracted from the model output 
# Figures 2 and 3 used interaction estimates which have been rescaled with a pop dyn model

# get real data
fecundities <- read.csv('2.case_study/data/fecundities0.csv', stringsAsFactors = F)
seeds <- fecundities$seeds
focalobs <- as.numeric(as.factor(fecundities$focal))

#################################################################################
# MAIN TEXT Figure 1: Posterior predictive check on the joint model 
##################################################################################

mu2 <- read.csv('2.case_study/model/output/mu2_samples.csv')
phi <- read.csv('2.case_study/model/output/disp_dev_samples.csv')
phi <- (phi^2)^(-1)

# select 1000 samples from the 80% posterior interval
mu2 <- apply(mu2, 2, function(p) {sample(p[p > quantile(p, 0.1) & p < quantile(p, 0.9)], size = 1000)})
phi <- apply(phi, 2, function(p) {sample(p[p > quantile(p, 0.1) & p < quantile(p, 0.9)], size = 1000)})
  
# generating posterior predictions
seed_pred <- matrix(nrow = dim(mu2)[1], ncol = dim(mu2)[2])
for (i in 1:dim(mu2)[1]) {     # for each posterior draw
  for (j in 1:dim(mu2)[2]) {    # for each observation 
    # draw from the predicted distribution
    seed_pred[i, j] <- rnbinom(1, mu = mu2[i, j], size = phi[i, focalobs[j]])  
  }
}
# log transform seed predictions
seed_pred <- log(seed_pred)
# and using just the median of the parameters
m.mu2 <- apply(mu2, 2, median)
m.phi <- apply(phi, 2, median)
mean_seed_pred <- sapply(1:length(m.mu2) , function(x) {
  rnbinom(1, mu = m.mu2[x], size = m.phi[focalobs[x]])
})
mean_seed_pred <- log(mean_seed_pred)

# get maximum density for plot limits
max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                     max(density(seeds)$y)))


png('3.analyses/figures_mss/postpredch_mu2.png', width=1563, height=1563)
par(cex = 4.5)
# start a plot with the first draw 
ppc.plot <- plot(density(seed_pred[1, ]), 
                 type = 'n',
                 bty = 'n',
                 ylim = c(0, max.density+0.03), 
                 xlim = c(-1, 10),
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

rm(mu2, ppc.plot, seed_pred)

##################################################################################
# MAIN TEXT Figure 2: Networks (Facilitative vs Competitive) 
##################################################################################

load('2.case_study/model/transformed/scaled_betaij_matrices.Rdata')

# Median interaction estimates for joint model
bij.med <- apply(scaled_betas, c(1, 2), median)
bij.med <- bij.med[ , 1:nrow(bij.med)]
bij.med <- t(bij.med)  # this is so qgraph points the arrows towards i

# set up colours for nodes (linking to next figure)
invasives <- c('ARCA', 'PEAI', 'HYPO')
foundation <- c('VERO', 'POCA')
keyst <- c('GITE')
all.sp <- rep('white', length(dimnames(bij.med)$species))
names(all.sp) <- dimnames(bij.med)$species
# all.sp[invasives] <- 'tomato'
all.sp[foundation] <- 'purple'
all.sp[keyst] <- 'aquamarine3'
all.sp[invasives] <- 'firebrick3'
# node.outline.width <- rep(1, length(dimnames(bij.med)$species))
# names(node.outline.width) <- dimnames(bij.med)$species
# node.outline.width[ c(invasives, foundation, keyst)] <- 3

library(qgraph)

png('3.analyses/figures_mss/networks_C_F.png', 
    width = 1875, height = 3750, units = 'px')
par(mfrow=c(2, 1))
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

dev.off()

##################################################################################
# MAIN TEXT Figure 3: Applications - species effects 
##################################################################################

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

# get the effects of species on others - FOCALS ONLY
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
keyst <- c('GITE')

# plot! 
png('3.analyses/figures_mss/species_effects.png', 
    width = 1563, height = 4688, units = 'px')
par(mfrow=c(3,1), cex=4)

# 1. intra vs abundance
#-------------------
plot(intras, sp.abunds.rep,
     xlab = 'Intraspecific interactions (self-regulation)',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2)
title(main = 'A', adj = 0)
polygon(x = c(-0.6 , 0, 0, -0.6), y = c(4.4, 4.4, 9.3, 9.3),
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
text(mean.intra[foundation], sp.abunds[foundation], 
     labels = names(mean.intra[foundation]), 
     pos = 4, col = 'darkorchid4', offset = 1, cex = 1.4)


# 2. out-strength vs abundance 
#-----------------------------
plot(sum.out, sp.abunds.rep, 
     xlab = 'Net effect on neighbours (Out-strength)',
     ylab = 'Log abundance', las = 1, type = 'n', bty = 'n', cex.lab = 1.2)
title(main = 'B', adj = 0)
polygon(x = c(-4 , 0, 0, -4), y = c(4.4, 4.4, 9.3, 9.3),
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
points(mean.sum.out[keyst], sp.abunds[keyst], pch = 24, 
       bg = 'chartreuse3', cex = 1.8)
lab.x.pos <- mean.sum.out[keyst]
text(lab.x.pos + .3, sp.abunds[keyst], 
     labels = names(mean.sum.out[keyst]), pos = 3, col = 'springgreen4', 
     cex = 1.4)


# 3. competitive vs facilitative effects
#-------------------------------------
plot(sum.comp, -sum.faci,  
     las = 1, bty = 'n', cex.lab = 1.2,
     # type = 'n',
     pch = 16, col = 'grey',
     xlab = 'Absolute sum of competitive effects', 
     ylab = 'Absolute sum of facilitative effects')
title(main = 'C', adj = 0)
abline(0, 1, cex=0.5, col = 'snow4')
abline(v=median(sum.comp), lty = 2)
abline(h=median(-sum.faci), lty = 2)
points(med.comp, -med.faci, 
       pch = 23, cex=1, bg='black')
points(med.comp[invasives], -med.faci[invasives],  
       pch = 24, bg = 'red', cex = 1.8)
lab.x.pos <- med.comp[invasives] 
lab.x.pos['PEAI'] <- lab.x.pos['PEAI'] - 0.15
lab.y.pos <- (-med.faci[invasives] + 0.11)
lab.y.pos['PEAI'] <- lab.y.pos['PEAI'] +0.05
lab.y.pos['ARCA'] <- lab.y.pos['ARCA'] - 0.2
text(lab.x.pos, lab.y.pos, 
     labels = names(med.comp[invasives]), 
     pos =  4, col = 'darkred', offset = 0.3, cex = 1.4)

dev.off()

##################################################################################
# SUPPS Figure 1: analyses of 1 gamma parameter
##################################################################################

library(rstan)

load('2.case_study/model/output/model_fit.Rdata')
sfit <- summary(fit, pars = 'gamma_i')$summary
names(which(sfit[ , 'Rhat'] == max(sfit[ , 'Rhat'])))

png('3.analyses/figures_mss/gamma_21_joint.png', width=800, height=500)
stan_par(fit, names(which(sfit[ , 'Rhat'] == max(sfit[ , 'Rhat']))))
dev.off()

rm(fit)

##################################################################################
# SUPPS Figure 2: Overlapping posterior density plots of gamma
##################################################################################

library(rethinking)
library(tidyverse)
library(ggjoy)
library(magrittr)

# extract gamma samples 
load('2.case_study/model/output/model_fit.Rdata')
joint_gamma <- extract.samples(fit, pars = 'gamma_i')[[1]]
rm(fit)

load('mc_nddm/model_fit.Rdata')
ndd_gamma <- extract.samples(fit, pars = 'gamma_i')[[1]]
rm(fit)

load('mc_rim/adapt_delta_099_treedepth_20/model_fit.Rdata')
rim_gamma <- extract.samples(fit, pars = 'gamma_i')[[1]]
rm(fit)

# set up plot! 
pl <- vector("list", length = 22)

for (i in 1:22) {
  
  Samples <- c(joint_gamma[ , i], ndd_gamma[ , i], rim_gamma[ , i])
  Model <-  c(rep('Joint', 8000), rep('NDDM', 8000), rep('RIM', 8000))
  df <- list(Model, Samples)
  df <- as.data.frame(df, col.names = c('Model', 'Samples'))
  
  df2 <- group_by(df, Model) %>% summarise('Samples' = median(Samples))
  
  g <- ggplot(df, aes(x=Samples, y=Model))+
    geom_joy(scale = 2, alpha=0.5) +
    scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(name = NULL, expand=c(0.01, 0)) +
    geom_point(data = df2, col = 'red') +
    theme_joy()
  
  pl[[i]] <- g
  
}

png('3.analyses/figures_mss/gamma_i_all_models.png', width=800, height=1200)
cowplot::plot_grid(plotlist = pl, nrow = 6)
dev.off()


##################################################################################
# SUPPS Figure 3: Posterior predictive check on the RIM only 
##################################################################################

# Do the same as for Fig 1 but for mu

# extract mu 
mu <- read.csv('2.case_study/model/output/mu_samples.csv')
mu <- apply(mu, 2, function(p) {sample(p[p > quantile(p, 0.1) & p < quantile(p, 0.9)], size = 1000)})

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
# and using just the median of the parameters
m.mu <- apply(mu, 2, median)
m.phi <- apply(phi, 2, median)
mean_seed_pred <- sapply(1:length(m.mu) , function(x) {
  rnbinom(1, mu = m.mu[x], size = m.phi[focalobs[x]])
})
mean_seed_pred <- log(mean_seed_pred)

# get maximum density for plot limits
max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                     max(density(seeds)$y)))


png('3.analyses/figures_mss/postpredch_mu.png', width=500, height=500)
# start a plot with the first draw 
ppc.plot <- plot(density(seed_pred[1, ]), 
                 type = 'n',
                 bty = 'n',
                 ylim = c(0, max.density+0.03),  
                 xlim = c(-1, 10),
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

##################################################################################
# SUPPS Figure 4: IFM vs RIM estimates 
##################################################################################

# Load data
rim_mat <- as.matrix(read.csv('2.case_study/model/output/RIM_betaij_samples.csv'))
nddm_mat <- as.matrix(read.csv('2.case_study/model/output/NDDM_betaij_samples.csv'))
# sample from 80% confidence interval
rim_mat <- apply(rim_mat, 2, function(p) {sample(p[p > quantile(p, 0.1) & 
                                                           p < quantile(p, 0.9)], size = 1000)})
nddm_mat <- apply(nddm_mat, 2, function(p) {
  if (sum(p) != 0) {sample(p[p > quantile(p, 0.1) & p < quantile(p, 0.9)], size = 1000)} else {rep(0, 1000)}})

mn <- min(nddm_mat)
mx <- max(nddm_mat)

# plot!
png('3.analyses/figures_mss/interaction_estimates.png', width = 1400, height = 700)
par(mfrow=c(2,2), cex = 1.5)
hist(nddm_mat[ , apply(nddm_mat, 2, function(x) {all(x != 0)})], 
     xlab = "", breaks = 30,
     main = "Realised interactions (NDDM)", 
     xlim = c(mn, mx))
plot(rim_mat[nddm_mat!=0], nddm_mat[nddm_mat!=0], 
     ylab = 'NDDM estimates', xlab='Response*Impact estimates',
     xlim = c(mn, mx), ylim = c(mn, mx), 
     pch = 16, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
abline(0,1)
hist(rim_mat[ , apply(nddm_mat, 2, function(x) {all(x == 0)})],  
     xlab = "", breaks = 30,
     main = 'Unrealised interactions (RIM)', 
     xlim = c(mn, mx))
hist(rim_mat[ , apply(nddm_mat, 2, function(x) {all(x != 0)})],  
     xlab = "", breaks = 30,
     main = 'Realised interactions (RIM)', 
     xlim = c(mn, mx))
dev.off()


# # # Potential figure for Supps ?
# # #------------------------------
# response <- as.matrix(read.csv('2.case_study/model/output/response_samples.csv'))
# effect <- as.matrix(read.csv('2.case_study/model/output/effect_samples.csv'))
# S <- dim(response)[2]
# 
# png('3.analyses/figures_mss/response_impact.png', width = 400, height = 700)
# plot(response, effect[ , 1:S], type = 'n',
#      xlab = expression(italic('response i')), ylab = expression(italic('impact i')))
# abline(h = 0, lty = 2)
# points(response, effect[ , 1:S],
#      pch = 16,
#      col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2))
# points(apply(response, 2, median), apply(effect[ , 1:S], 2, median),
#        pch = 18, cex = 2)
# dev.off()
