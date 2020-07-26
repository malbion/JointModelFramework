setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')

load('../2018_Compnet/stormland/model/output/0/post_draws.Rdata')

fit_sum <- read.csv('3.case_study/model/output/summary_of_draws.csv', row.names = 1)
key_speciesID <- unlist(read.csv('3.case_study/data/key_speciesID0.csv', stringsAsFactors = F))

# Interactions (alphas)
# get posterior draws for all interactions (from the IFM)
alphas <- joint.post.draws$ifm_alpha
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

rm(joint.post.draws)

png(paste0('model/validation/', comm, '/ifm_vs_rem_alphas.png'))
plot(alphas[alphas!=0], re_alphas[alphas!=0], xlab = 'IFM betas', ylab='Response*Effect estimates',
     xlim = c(min(re_alphas), max(re_alphas)),
     ylim = c(min(re_alphas), max(re_alphas)), pch = 16, 
     col = rgb(red = 1, green = 1, blue = 1, alpha = 0.4))
abline(0,1)
dev.off()

unobs <- re_alphas[ , apply(alphas, 2, function(x) {all(x == 0)})]

png('2.analyses/figures/alpha_est_distr.png', width = 500, height = 500)
par(mfrow=c(3,1))
hist(alphas[ , apply(alphas, 2, function(x) {all(x != 0)})], xlab = "", breaks = 30,
     main = "IFM interaction estimates (0's removed)", xlim = c(min(re_alphas), max(re_alphas)))
hist(re_alphas[ , apply(alphas, 2, function(x) {all(x != 0)})],  xlab = "", breaks = 30,
     main = 'REM interaction estimates - Unobserved interactions removed', xlim = c(min(re_alphas), max(re_alphas)))
hist(unobs,  xlab = "", main = 'REM interaction estimates - Unobserved only', breaks = 30,
     xlim = c(min(re_alphas), max(re_alphas)))
dev.off()

png('2.analyses/figures/interaction_ests.png', width = 1400, height = 700)
par(mfrow=c(2,2), cex = 1.5)
hist(alphas[ , apply(alphas, 2, function(x) {all(x != 0)})], xlab = "", breaks = 30,
     main = "IFM interaction estimates (0's removed)", xlim = c(min(re_alphas), max(re_alphas)))
plot(re_alphas[alphas!=0], alphas[alphas!=0], ylab = 'IFM betas', xlab='Response*Effect estimates',
     xlim = c(min(re_alphas), max(re_alphas)),
     ylim = c(min(re_alphas), max(re_alphas)), pch = 16, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
abline(0,1)
hist(unobs,  xlab = "", main = 'REM interaction estimates - Unobserved only', breaks = 30,
     xlim = c(min(re_alphas), max(re_alphas)))
hist(re_alphas[ , apply(alphas, 2, function(x) {all(x != 0)})],  xlab = "", breaks = 30,
     main = 'REM interaction estimates - Unobserved interactions removed', xlim = c(min(re_alphas), max(re_alphas)))
dev.off()

# Find estimates from REM / IFM
#-------------------------------
load('3.case_study/model/transformed/scaled_alpha_matrices.Rdata')
alpha_mat <- apply(scaled_alphas, c(1, 2), mean)
rm(scaled_alphas)

obs <- read.csv('3.case_study/data/obs_interact0.csv')
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

png('2.analyses/figures/ifm_only.png', width = 600, height = 240, units = 'px')
qgraph(ifm_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
dev.off()

png('2.analyses/figures/rem_only.png', width = 600, height = 240, units = 'px')
qgraph(rem_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
dev.off()

png('2.analyses/figures/all_vs_rem.png', width = 1400, height = 1040, units = 'px')
par(mfrow=c(2,1))
qgraph(alpha_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle',
       title = 'A', title.cex =5)
qgraph(rem_mat, negCol = 'royalblue4', posCol = 'orange', layout = 'circle',
       title = 'B', title.cex =5)
dev.off()

