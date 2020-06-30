setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')

load('3.case_study/model/transformed/scaled_alpha_matrices.Rdata')
# Calculate metrics
source('../2018_Compnet/stormland/functions/analyses/calc_sp_ntwk_metrics.R')
sp.abunds <- read.csv('3.case_study/data/plot_species_abundances.csv', stringsAsFactors = F)
metrics <- calc.sp.ntwk.metrics(0, scaled_alphas, sp.abunds)

write.csv(metrics, '2.analyses/species_metrics_mat.csv', row.names = F)
metrics <- read.csv('2.analyses/species_metrics_mat.csv')
met.sum <- apply(metrics, c(1,2), mean)


source('../2018_Compnet/stormland/functions/analyses/calc_com_ntwk_metrics.R')


