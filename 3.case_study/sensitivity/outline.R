# Checking model sensitivity

# 13.6% of pairwise interactions were unobserved 
# run model with 20%, 30%, 40% and 50% unobserved interactions 
# 100 sims per category, interactions randomly removed


# plot estimates for the removed interactions vs their 'true' value for each category
# plot the sd of alpha estimates vs the number of observations for that interaction


# Process for removing observations: 
# - randomly select focal observations to remove in blocks of 10 observations
# - or select an interaction to remove and remove all focal observations corresponding to that interactions? 
# - calculate resulting number of observations for each interaction
# - and calculate resulting number of missing observations 


setwd('Dropbox/Work/Projects/2020_Methods_for_compnet/3.case_study/sensitivity/')


fecundities <- read.csv('../data/fecundities0.csv', stringsAsFactors = F)
# create a matrix which tallies the number of observations for each focal and neighbour
fec <- fecundities[ , -c(1:4)]
fec[fec>0] <- 1
fec <- split(fec, fecundities$focal)
obs <- do.call(rbind, lapply(fec, colSums))
library(reshape2)
obs_long <- melt(obs)


# let's randomly remove alpha_ARCA_EROD

arca_fec <- fecundities[fecundities$focal == 'ARCA' & fecundities$EROD != 0, ]
rownames(arca_fec)

fec2 <- fecundities[-as.numeric(rownames(arca_fec)), ]

# urgh... alternatively just create some fake data?


