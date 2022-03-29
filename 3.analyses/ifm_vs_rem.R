# Plot networks of NDDM and RIM interaction estimates
#------------------------------------------------------

library(here)
setwd(here())

keyspeciesID <- unlist(read.csv('2.case_study/data/key_speciesID0.csv', stringsAsFactors = F))
keyneighbourID <- unlist(read.csv('2.case_study/data/key_neighbourID0.csv', stringsAsFactors = F))

# Get raw (unscaled) interaction estimates
joint_betaij <- read.csv('2.case_study/model/output/joint_betaij_samples.csv')
rim_betaij <- read.csv('2.case_study/model/output/RIM_betaij_samples.csv')
nddm_betaij <- read.csv('2.case_study/model/output/NDDM_betaij_samples.csv')

# get median species x species interaction matrices 
joint_mat <- matrix(apply(joint_betaij, 2, median), byrow = T,
                    nrow = length(keyspeciesID), ncol = length(keyneighbourID))
colnames(joint_mat) <- keyneighbourID
rownames(joint_mat) <- keyspeciesID
rim_mat <- matrix(apply(rim_betaij, 2, median), byrow = T,
                    nrow = length(keyspeciesID), ncol = length(keyneighbourID))
colnames(rim_mat) <- keyneighbourID
rownames(rim_mat) <- keyspeciesID
nddm_mat <- matrix(apply(nddm_betaij, 2, median), byrow = T,
                    nrow = length(keyspeciesID), ncol = length(keyneighbourID))
colnames(nddm_mat) <- keyneighbourID
rownames(nddm_mat) <- keyspeciesID


library(qgraph)
qgraph(joint_mat[ , 1:length(keyspeciesID)], negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
qgraph(nddm_mat[ , 1:length(keyspeciesID)], negCol = 'royalblue4', posCol = 'orange', layout = 'circle')
qgraph(rim_mat[ , 1:length(keyspeciesID)], negCol = 'royalblue4', posCol = 'orange', layout = 'circle')

png('3.analyses/figures/joint_network.png', width = 1200, height = 480, units = 'px')
qgraph(t(joint_mat[ , 1:length(keyspeciesID)]), negCol = 'royalblue4', posCol = 'orange', 
       layout = 'circle', maximum = 1.6, diag = T)
dev.off()

png('3.analyses/figures/NDDM_network.png', width = 1200, height = 480, units = 'px')
qgraph(t(nddm_mat[ , 1:length(keyspeciesID)]), negCol = 'royalblue4', posCol = 'orange', layout = 'circle', maximum = 1.6)
dev.off()

png('3.analyses/figures/RIM_network.png', width = 1200, height = 480, units = 'px')
qgraph(t(rim_mat[ , 1:length(keyspeciesID)]), negCol = 'royalblue4', posCol = 'orange', layout = 'circle', maximum = 1.6)
dev.off()

rim_mat_infer <- rim_mat
rim_mat_infer[which(nddm_mat == 0)] <- 0  # only keep the interactions which are also inferred by the NDDM
png('3.analyses/figures/RIM_inferrables_network.png', width = 1200, height = 480, units = 'px')
qgraph(t(rim_mat_infer[ , 1:length(keyspeciesID)]), negCol = 'royalblue4', posCol = 'orange', layout = 'circle', maximum = 1.6)
dev.off()

rim_mat_noinfer <- rim_mat
rim_mat_noinfer[which(nddm_mat != 0)] <- 0 # only keep the interactions which aren't inferred by the NDDM
png('3.analyses/figures/RIM_noninferrables_network.png', width = 1200, height = 480, units = 'px')
qgraph(t(rim_mat_noinfer[ , 1:length(keyspeciesID)]), negCol = 'royalblue4', posCol = 'orange', layout = 'circle', maximum = 1.6)
dev.off()

