

# Case study for Methods - on Trace's 'full' community dataset

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')

# How many REM estimates did we use? 
#-----------------------------------
# for focals x focals: 
obs <- read.csv('3.case_study/data/obs_interact0.csv', stringsAsFactors = F)
table(obs$cooc_code)
410/484*100
# for the whole thing: check stan.data$I ! 
687/(22*52)*100

# Alpha stuff
load('3.case_study/model/transformed/scaled_alpha_matrices.Rdata')
range(scaled_alphas)
# Do stuff on the alphas
length(scaled_alphas[scaled_alphas > 0])/length(scaled_alphas)*100
length(scaled_alphas[scaled_alphas < 0])/length(scaled_alphas)*100

# number of neighbour a focal competes with vs helps? 
# total N they have a competitive effect on 
compef <- apply(scaled_alphas, 3, function(sm){
  apply(sm[ , 1:22], 2, function(r){length(r[r>0])})
})
compef <- t(compef)
colnames(compef) <- dimnames(scaled_alphas)[[1]]
colMeans(compef)

# total N they have a facilitative effect on 
facef <- apply(scaled_alphas, 3, function(sm){
  apply(sm[ , 1:22], 2, function(r){length(r[r<0])})
})
facef <- t(facef)
colnames(facef) <- dimnames(scaled_alphas)[[1]]
colMeans(facef)


# get species metrics
source('2.analyses/calc_sp_ntwk_metrics.R')
sp.abunds <- read.csv('3.case_study/data/plot_species_abundances.csv', stringsAsFactors = F)
smm <- calc.sp.ntwk.metrics('0', scaled_alphas, sp.abunds)
require(plyr)
smm <- adply(smm, c(1, 3))

smmfoc <- smm[is.na(smm$aii) == F, ]
smmfoc$species <- droplevels(smmfoc$species)

lsf <- split(smmfoc, smmfoc$species)
parofinterest <- c('sum_aij', 'sum_aji', 'C_sum_aij', 'C_sum_aji', 'F_sum_aij', 'F_sum_aji')
species.chars <- lapply(lsf, function(sp){
  means <- colMeans(sp[ , 3:23])
  sds <- apply(sp[ , 3:23], 2, sd)
  names(sds) <- paste0('sd_', names(sds))
  temp <- c(unique(sp$species), unique(sp$com.abund), unique(sp$rel.abund), means, sds)
  names(temp)[1:3] <- c('species', 'com.abund', 'rel.abund')
  return(temp)
  })
species.chars <- do.call(rbind, species.chars) 
species.chars <- cbind(species.chars, colMeans(compef))
species.chars <- cbind(species.chars, colMeans(facef))
colnames(species.chars)[46:47] <- c('N_compef', 'N_facef')

                                    