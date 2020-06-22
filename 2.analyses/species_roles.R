setwd('~/Dropbox/Work/Projects/2018_Compnet/stormland/')

load('model/transformed/0/scaled_alpha_matrices.Rdata')
# Calculate metrics
source('functions/analyses/calc_sp_ntwk_metrics.R')
sp.abunds <- read.csv('clean_data/plot_species_abundances.csv', stringsAsFactors = F)
metrics <- calc.sp.ntwk.metrics(0, scaled_alphas, sp.abunds)

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/2.analyses/')
write.csv(metrics, 'species_metrics_mat.csv', row.names = F)

met.sum <- apply(metrics, c(1,2), mean)


source('functions/analyses/calc_com_ntwk_metrics.R')


