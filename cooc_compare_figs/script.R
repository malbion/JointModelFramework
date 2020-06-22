setwd('~/Dropbox/Work/Projects/2018_Compnet/stormland/')
require(qgraph)


inter.mat <- list()
inter.mat <- sapply(c(1, 2, 3), function(c){
load(paste0('model/transformed/', c, '/scaled_alpha_matrices.Rdata'))
inter.mat[[as.numeric(c)]] <- scaled_alphas
})


inter <- lapply(inter.mat, function(x) {apply(x, c(1, 2), mean)})
rm(inter.mat)
library(qgraph)


cooc <- sapply(c(1, 2, 3), function(com) {
read.csv(paste0('../3.results/archive/2.co_occurrence/', ccc[com], '_cooc.csv'), stringsAsFactors = F, row.names = 1)
}, simplify = F, USE.NAMES =T)


sapply(c(1, 2, 3), function(x){
  
  png(paste0('../../2020_Methods_for_compnet/cooc_compare_figs/subntwk', x, '.png'), 
      width = 900, height = 720, units = 'px')
  par(mfrow = c(2, 1))
  qgraph((-inter[[x]][ , 1:nrow(inter[[x]])]))
  qgraph(cooc[[x]][rownames(inter[[x]]), colnames(cooc[[x]]) %in% rownames(inter[[x]])])
  dev.off()
})