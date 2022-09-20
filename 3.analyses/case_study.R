# Quantitative analysis of the case study interaction matrices

library(here)
setwd(here())


load('2.case_study/model/transformed/scaled_betaij_matrices.Rdata')

range(scaled_betas) # range of all strengths 

# CALCULATION ON WHOLE B MATRIX 
#------------------------------
# # BELOW INCLUDES INTRASPECIFIC INTERACTIONS IN THE CALCULATIONS
# # percentage of competitive interactions
# length(scaled_betas[scaled_betas > 0])/length(scaled_betas)*100
# # percentage of facilitative interactions
# length(scaled_betas[scaled_betas < 0])/length(scaled_betas)*100
# # interactions between focal species only
# foc_only <- apply(scaled_betas, 3, function(smp){smp <- smp[1:22, 1:22]})
# # percentage of focal-focal competitive interactions
# length(foc_only[foc_only > 0])/length(foc_only)*100
# # percentage of focal-focal facilitative interactions
# length(foc_only[foc_only < 0])/length(foc_only)*100

# # BELOW DOES NOT INCLUDE INTRASPECIFIC INTERACTIONS
# # FOCAL-NEIGHBOUR interactions with NO INTRAS
# focneigh <- apply(scaled_betas, 3, function(mat) {diag(mat) <- 0; return(mat)})
# focneigh <- focneigh[focneigh != 0]
# # percentage of competitive interactions
# length(focneigh[focneigh > 0])/length(focneigh)*100
# # percentage of facilitative interactions
# length(focneigh[focneigh < 0])/length(focneigh)*100
# # FOCAL-FOCAL interactions with NO INTRAS
# focfoc <- apply(scaled_betas, 3, function(mat) {diag(mat) <- 0; return(mat[1:22, 1:22])})
# focfoc <- focfoc[focfoc != 0]
# # percentage of competitive interactions
# length(focfoc[focfoc > 0])/length(focfoc)*100
# # percentage of facilitative interactions
# length(focfoc[focfoc < 0])/length(focfoc)*100

# CALCULATIONS ON MEDIAN INTERACTION STRENGTHS
#---------------------------------------------
med_betas <- apply(scaled_betas, c(1,2), median)
# Get intras
med_bii <- diag(med_betas)
# remove intras from calculations
diag(med_betas) <- 0
# get focal x neighbour matrix
med_betas_FN <- med_betas[med_betas!=0]
# percentage of competitive interactions
length(med_betas_FN[med_betas_FN>0])/length(med_betas_FN)*100
# percentage of facilitative interactions
length(med_betas_FN[med_betas_FN<0])/length(med_betas_FN)*100
# get focal x focal matrix
med_betas_FF <- med_betas[1:22, 1:22]
med_betas_FF <- med_betas_FF[med_betas_FF!=0]
# percentage of competitive interactions
length(med_betas_FF[med_betas_FF>0])/length(med_betas_FF)*100
# percentage of facilitative interactions
length(med_betas_FF[med_betas_FF<0])/length(med_betas_FF)*100
# look at intras only 
length(med_bii[med_bii>0])/length(med_bii)*100  # competitive
length(med_bii[med_bii<0])/length(med_bii)*100  # facilitative

# Asymmetry 
med_betas_FF <- med_betas[1:22, 1:22]
asym <- matrix(data = NA, nrow = 22, ncol = 22)
for (i in 1:nrow(med_betas_FF)) {
  for (j in 1: ncol(med_betas_FF)) {
    # return 1 if asymmetric, 0 if they match sign
    asym[i , j] <- ifelse(sign(med_betas_FF[i , j]) == sign(med_betas_FF[j, i]), 0, 1)
  }
}
asym <- asym[upper.tri(asym)]
length(asym[asym == 1])/length(asym)*100

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
