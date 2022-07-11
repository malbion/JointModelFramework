# Quantitative analysis of the case study interaction matrices

library(here)
setwd(here())


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
write.csv(smm, '3.analyses/smm.csv', row.names = F)

smmii <- smm[is.na(smm$aii) == F, ]
median(smmii$perc.exploiter)*2

library(magrittr)
library(dplyr)

# intraspecific interactions - overlap with 0
smmii$species <- droplevels(smmii$species)
smmii %>% group_by(species) %>% summarise(min = min(aii),
                                          max = max(aii)) -> foo
foo[which(sign(foo$min) == sign(foo$max)), ]

# interspecific interactions - overlap with 0
minval <- apply(scaled_betas, c(1,2), min)
maxval <- apply(scaled_betas, c(1,2), max)
# remove intraspecific interactions 
diag(minval) <- NA
diag(maxval) <- NA
# proportion which do not overlap with 0
length(minval[which(sign(minval) == sign(maxval))]) / (dim(minval)[1] * dim(minval)[2] - 22) * 100
# check on focal x focal only 
minval <- minval[ , 1:nrow(minval)]
maxval <- maxval[ , 1:nrow(maxval)]
