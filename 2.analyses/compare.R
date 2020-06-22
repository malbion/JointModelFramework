

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/2.analyses/')

load('../../2018_Compnet/stormland/model/transformed/0/scaled_alpha_matrices.Rdata')

alpha_mat <- apply(scaled_alphas, c(1, 2), mean)

#-----------------------------------
# Compare to co-occurrence network |
#-----------------------------------
require(cooccur)
require(reshape2)
# 1. Build a cooccurrence network
data_folder <- '~/Dropbox/Work/Projects/2018_Compnet/1.data/clean/'
input_folder <- '~/Dropbox/Work/Projects/2018_Compnet/3.results/1.joint_model'

fecundities <- read.csv(paste0('../../2018_Compnet/stormland/clean_data/fecundities', comm, '.csv'), stringsAsFactors = F)
key_speciesID <- unlist(read.csv(paste0('../../2018_Compnet/stormland/clean_data/key_speciesID', comm, '.csv'), stringsAsFactors = F))
key_neighbourID <- unlist(read.csv(paste0('../../2018_Compnet/stormland/clean_data/key_neighbourID', comm, '.csv'), stringsAsFactors = F))

# each observation is a neighbourhood, hence a 'plot'
for(i in unique(fecundities$focal)) {
  fecundities[ , i] <- fecundities[ , i] + as.numeric(fecundities$focal == i)
} # add focal observation to the neighbourhood abundances
# we now have a species.site matrix
spp.site <- t(fecundities[ , -c(1:4)])

# this function applies the probabilistic model of species cooccurrence (Veech 2013)
cooc.site.sp <- cooccur(spp.site, spp_names = T,  
                        thresh = F, #if TRUE cooc < 1 are removed
)

# 2. Extract the community matrix
associations <- effect.sizes(cooc.site.sp)
assocs <- dcast(associations, sp1 ~ sp2)
assocs <- as.data.frame(assocs[ , -1], row.names = as.character(assocs[ , 1]))
# Why is this missing a row and column? 
# This seems to be the same row and column across both 0to15 and 15plus...
# row = others, column = CHEI
assocs <- assocs[c(key_speciesID, sort(key_neighbourID[!key_neighbourID %in% key_speciesID])), ]
# convert 'NA' rowname to 'others
rownames(assocs) <- c(key_speciesID, sort(key_neighbourID[!key_neighbourID %in% key_speciesID]))
# add a column for the missing species
assocs[ , key_neighbourID[!key_neighbourID %in% colnames(assocs)]] <- NA
assocs <- assocs[, c(key_speciesID, sort(key_neighbourID[!key_neighbourID %in% key_speciesID]))]
s <- assocs
s[lower.tri(s)] = t(s)[lower.tri(s)]
p <- assocs
p[upper.tri(p)] = t(p)[upper.tri(p)]
s[is.na(s)] <- 0
p[is.na(p)] <- 0
assocs <- s + p
assocs <- as.matrix(assocs)
# return intra-sp associations as NAs (cooccur can't estimate thos)
diag(assocs) <- NA

write.csv(assocs, paste0(comm, '_cooc.csv'))

png(paste0('../../2020_Methods_for_compnet/cooc_compare_figs/subntwk', x, '.png'), 
    width = 900, height = 720, units = 'px')
par(mfrow = c(2, 1))
qgraph((-alpha_mat[ , 1:nrow(alpha_mat)]))
qgraph(assocs[rownames(alpha_mat), colnames(assocs) %in% rownames(alpha_mat)])
dev.off()
