# Plot networks of NDDM and RIM interaction estimates
#------------------------------------------------------

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')

keyspeciesID <- unlist(read.csv('2.case_study/data/key_speciesID0.csv', stringsAsFactors = F))
keyspeciesID <- unlist(read.csv('2.case_study/data/key_neighbourID0.csv', stringsAsFactors = F))

# Get raw (unscaled) interaction estimates
joint_betaij <- read.csv('2.case_study/model/output/joint_betaij_samples.csv')
rim_betaij <- read.csv('2.case_study/model/output/RIM_betaij_samples.csv')
nddm_betaij <- read.csv('2.case_study/model/output/NDDM_betaij_samples.csv')

# get median species x species interaction matrices 


# load('2.case_study/model/output/post_draws.Rdata')
# 
# fit_sum <- read.csv('2.case_study/model/output/summary_of_draws.csv', row.names = 1)
# key_speciesID <- unlist(read.csv('2.case_study/data/key_speciesID0.csv', stringsAsFactors = F))

# Interactions (alphas)
# get posterior draws for all interactions (from the IFM)


# Find estimates from REM / IFM
#-------------------------------
load('2.case_study/model/transformed/scaled_betaij_matrices.Rdata')
alpha_mat <- apply(scaled_alphas, c(1, 2), mean)
rm(scaled_alphas)

obs <- read.csv('2.case_study/data/obs_interact0.csv')
obs <- obs[obs$neighb_is_focal == T, ]
library(reshape2)
obs <- dcast(obs, focal ~ neighbour, value.var = 'total_density')
rownames(obs) <- obs[ , 1]
obs <- obs[ , 2:23]

alpha_mat <- alpha_mat[ , 1:nrow(alpha_mat)]
ifm_mat <- alpha_mat
rem_mat <- alpha_mat

ifm_mat[which(obs==0)] <- 0
rem_mat[which(obs>0)] <- 0

library(qgraph)
qgraph(alpha_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
qgraph(ifm_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
qgraph(rem_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')

png('3.analyses/figures/ifm_only.png', width = 600, height = 240, units = 'px')
qgraph(ifm_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
dev.off()

png('3.analyses/figures/rem_only.png', width = 600, height = 240, units = 'px')
qgraph(rem_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
dev.off()

png('3.analyses/figures/all_vs_rem.png', width = 1400, height = 1040, units = 'px')
par(mfrow=c(2,1))
qgraph(alpha_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle',
       title = 'A', title.cex =5)
qgraph(rem_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle',
       title = 'B', title.cex =5)
dev.off()

