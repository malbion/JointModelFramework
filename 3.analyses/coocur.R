

setwd('~/Dropbox/Work/Projects/2020_Methods_for_compnet/')

load('3.case_study/model/transformed/scaled_betaij_matrices.Rdata')

betaij_mat <- apply(scaled_betas, c(1, 2), median)

#-----------------------------------
# Compare to co-occurrence network |
#-----------------------------------
require(cooccur)
require(reshape2)
# 1. Build a cooccurrence network
data_folder <- '~/Dropbox/Work/Projects/2018_Compnet/1.data/clean/'
input_folder <- '~/Dropbox/Work/Projects/2018_Compnet/3.results/1.joint_model'
comm <- 0
fecundities <- read.csv(paste0('3.case_study/data/fecundities', comm, '.csv'), stringsAsFactors = F)
key_speciesID <- unlist(read.csv(paste0('3.case_study/data/key_speciesID', comm, '.csv'), stringsAsFactors = F))
key_neighbourID <- unlist(read.csv(paste0('3.case_study/data/key_neighbourID', comm, '.csv'), stringsAsFactors = F))

# each observation is a neighbourhood, hence a 'plot'
for(i in unique(fecundities$focal)) {
  fecundities[ , i] <- fecundities[ , i] + as.numeric(fecundities$focal == i)
} # add focal observation to the neighbourhood abundances
# we now have a species.site matrix
spp.site <- t(fecundities[ , -c(1:4)])

# this function applies the probabilistic model of species cooccurrence (Veech 2013)
cooc.site.sp <- cooccur(spp.site, spp_names = T,  
                        thresh = F #if TRUE cooc < 1 are removed
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

write.csv(assocs, paste0('2.analyses/', comm, '_cooc.csv'))

library(qgraph)

# comp <- -alpha_mat[ , 1:nrow(alpha_mat)]
# comp[comp>0] <- 0
# facil <- -alpha_mat[ , 1:nrow(alpha_mat)]
# facil[facil<0] <- 0
# png('2.analyses/figures/comp_facil.png', width = 1400, height = 1040, units = 'px')
# par(mfrow = c(2, 1))
# qgraph(comp, posCol = 'royalblue4', negCol = 'orange', layout = 'circle',
#        title = 'A', title.cex =5)
# qgraph(facil, posCol = 'royalblue4', negCol = 'orange', layout = 'circle',
#        title = 'B', title.cex =5)
# dev.off()
# # qgraph((-alpha_mat[ , 1:nrow(alpha_mat)]))
# png('2.analyses/figures/coocur.png', width = 1100, height = 520, units = 'px')
# qgraph(assocs[rownames(alpha_mat), colnames(assocs) %in% rownames(alpha_mat)],
#        posCol = 'royalblue4', negCol = 'orangered', layout = 'circle',
#        title = 'C', title.cex =5)
# dev.off()
# 
# png('2.analyses/figures/all_inter.png', width = 600, height = 240, units = 'px')
# qgraph(-alpha_mat[ , 1:nrow(alpha_mat)], posCol = 'royalblue4', negCol = 'orange', layout = 'circle')
# dev.off()

smm <- read.csv('2.analyses/smm.csv', stringsAsFactors = F)
smmfoc <- smm[is.na(smm$aii) == F, ]
mean(smmfoc$perc.coop)
mean(smmfoc$perc.fullcomp)
mean(smmfoc$perc.exploited)
mean(smmfoc$perc.exploiter)


p2 = mean(smmfoc$perc.coop)
p1 = mean(smmfoc$perc.exploited) + mean(smmfoc$perc.exploiter)
p3 = mean(smmfoc$perc.fullcomp)

plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n",
     xlab="", ylab="", bty='n', las = 1)
polygon(c(0, 0.1, 0.1, 0), c(0, 0, p1, p1), col = 'orange', border = NA)
polygon(c(0, 0.1, 0.1, 0), c(p1, p1, p1+p2, p1+p2), col = 'grey90', border = NA)
polygon(c(0, 0.1, 0.1, 0), c(p1+p2, p1+p2, p1+p2+p3, p1+p2+p3), col = 'red', border = NA)
polygon(c(0, 0.1, 0.1, 0), c(0, 0, 1, 1))
text(c(0.15, 0.15, 0.15), y = c(0.8, 0.44, 0.2), c('-/-', '+/+', '-/+'), cex = 2)





