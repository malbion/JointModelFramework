# Quantitative analysis of the case study interaction matrices

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')


load('2.case_study/model/transformed/scaled_betaij_matrices.Rdata')

range(scaled_betas) # range of strengths 
# percentage of competitive interactions
length(scaled_betas[scaled_betas > 0])/length(scaled_betas)*100
# percentage of facilitative interactions
length(scaled_betas[scaled_betas < 0])/length(scaled_betas)*100
# interactions between focal species only
foc_only <- apply(scaled_betas, 3, function(smp){smp <- smp[1:22, 1:22]})
# percentage of focal-focal competitive interactions
length(foc_only[foc_only > 0])/length(foc_only)*100
# percentage of focal-focal facilitative interactions
length(foc_only[foc_only < 0])/length(foc_only)*100


# SUMMARISE SPECIES METRICS

# get species metrics
source('3.analyses/calc_sp_ntwk_metrics.R')
sp.abunds <- read.csv('2.case_study/data/plot_species_abundances.csv', stringsAsFactors = F)
smm <- calc.sp.ntwk.metrics('0', scaled_betas, sp.abunds)
require(plyr)
smm <- adply(smm, c(1, 3))
write.csv(smm, '2.analyses/smm.csv', row.names = F)

smm <- read.csv('2.analyses/smm.csv', stringsAsFactors = F)
                                    