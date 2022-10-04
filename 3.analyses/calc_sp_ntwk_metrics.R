# Malyon Bimler
# CompNet 2019 

# Species-level network analysis

# Calculating Effective Competition
#----------------------------------
require(reshape2)

efc <- function(imat, com, intras, ...) {
  fecund <- read.csv(paste0('2.case_study/data/fecundities', com, '.csv'), stringsAsFactors = F)
  fecund <- fecund[ , !colnames(fecund) %in% c('plot', 'seeds', 'focalID')]
  
  effective.comp <- sapply(dimnames(imat)[['species']], function(foc) {
    
    abunds <- fecund[fecund$focal == foc, ]
    abunds <- colMeans(abunds[, -1])
    alphs <-imat[foc, , ]
    
    # intraspecific
    ints <- mean(fecund[fecund$focal == foc, foc]) # mean intra abundance
    
    eff.comp <- rbind(abunds[rownames(alphs)]*alphs, ints*intras[foc, ])

    return(eff.comp)
  }, simplify = F, USE.NAMES = T)
  
  eff.c.aii <- lapply(effective.comp, function(x){x[ dim(x)[[1]], ]})
  eff.c.aii <- do.call(rbind, eff.c.aii)
  
  effective.comp <- lapply(effective.comp, function(x){colSums(x[ -dim(x)[[1]], ])})
  effective.comp <- do.call(rbind, effective.comp)
  
  res <- list(effective.comp, eff.c.aii)
  names(res) <- c('eff.comp', 'eff.comp.intra')
  return(res)
}

# 2. Given a species x neighbour x samples array of scaled interaction strengths (alpha) 
# return a df containing species-level network metrics 

calc.sp.ntwk.metrics <- function(
  com, # community
  imat, # species x neighbour x samples array of scaled interaction strengths
  abund, # plot-level abundances for species
                                 ...){
  
  # Get the number of species, neighbours and samples
  Sp <- as.numeric(dim(imat)[[1]])
  Nb <- as.numeric(dim(imat)[[2]])
  Sm <- as.numeric(dim(imat)[[3]])
  
  # check that neighbours are ordered with focals first 
  if(identical(dimnames(imat)[[1]], 
               dimnames(imat)[[2]][1:Sp]) == F){
    message('ERROR! Neighbours not ordered by focals first, results incorrect')
  }
  # set up results df 
  metricnames <- c('sum_aij', 'sum_|aij|', 'sum_aji', 'sum_|aji|',
                   'C_sum_aij', 'C_sum_aji', 'F_sum_aij', 'F_sum_aji', 
                   'aii', 'aii>0', # 'imbalance',  'indirect.effs', 
                   'H.in', 'H.out',  # 'rH.out', 'rH.in', 
                   'perc.exploited', 'perc.exploiter', 'perc.fullcomp', 
                   'perc.coop','comp.rank', 'response', 'effect', 
                   'eff.comp', 'eff.comp.intra', 'com.abund', 
                   'rel.abund', 'com')
  metrics <- array(data = NA, dim = c(Nb, length(metricnames), Sm),
                   dimnames = list('species' = dimnames(imat)[[2]], 
                                   'metrics' = metricnames,
                                   'samples' = dimnames(imat)[[3]]))
  
  # Extract intraspecific alphas from the community matrix 
  intras <- apply(imat, 3, diag)
  if(nrow(intras) != Sp){message('ERROR! Intraspecifics not properly extracted')}
  
  # Remove intras from the comm matrix for network metrics (they are treated differently)
  for(k in 1:Sm) {diag(imat[ , , k]) <- 0}
  
  # create a competition-only and a facilitation-only matrix to calculate the different metrics
  # on too 
  comp.only <- imat
  comp.only[comp.only < 0] <- 0
  facil.only <- imat
  facil.only[facil.only > 0] <- 0
  
  # STRENGTHS AND SENSITIVITIES
  #-----------------------------
  # in-strength
  metrics[ , 'sum_aij', ] <- rbind(apply(imat, 3, rowSums),
                                   matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  # in-sensitivity (aka in-strength on absolute alphas)
  metrics[ , 'sum_|aij|', ] <- rbind(apply(abs(imat), 3, rowSums),
                                     matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  # out-strength
  metrics[ , 'sum_aji', ] <- apply(imat, 3, colSums)
  # out-sensitivity (aka out-strength on absolute alphas)
  metrics[ , 'sum_|aji|', ] <- apply(abs(imat), 3, colSums)
  
  # same but on comp or facil only 
  metrics[ , 'C_sum_aij', ] <- rbind(apply(comp.only, 3, rowSums),
                                   matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'C_sum_aji', ] <- apply(comp.only, 3, colSums)
  metrics[ , 'F_sum_aij', ] <- rbind(apply(facil.only, 3, rowSums),
                                     matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'F_sum_aji', ] <- apply(facil.only, 3, colSums)
  
  # INTRASPECIFIC INTERACTIONS 
  metrics[ , 'aii', ] <- rbind(intras, matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  # competitive or facilitative? 
  metrics[ , 'aii>0', ] <- ifelse(metrics[ , 'aii', ] > 0, T, F) # 1 if competitive

  # Some of these are returning -Inf because C_sum_aij or aji is 0 - in these cases,
  # divide aii by 1 
  inf.to.rm <- which(metrics==-Inf | metrics==Inf, arr.ind = T)
  if(nrow(inf.to.rm) > 0) {
    for(n in 1:nrow(inf.to.rm)){
      metrics[inf.to.rm[n, 1],  inf.to.rm[n, 2], inf.to.rm[n, 3]] <- 
        metrics[inf.to.rm[n, 1], 'aii', inf.to.rm[n, 3]]
    }
  }

  # COMMUNITY-DEPENDANT METRICS
  #-----------------------------
  # # imbalance
  # metrics[ , 'imbalance', ] <- metrics[ , 'sum_aij', ] - metrics[ , 'sum_aji', ]
  
  
  # METRICS BETWEEN INTERACTING FOCALS ONLY (need both in and out strength)
  #---------------------------------------
  
  iemat <- sapply(1:Sm, function(n) {
    res <- matrix(data = NA, nrow = Sp, ncol = 6)
    
    for (i in 1:Sp) {
      # get index of neighbours which are also focals but remove species i 
      j <- 1:Sp
      j <- j[-i]
      
      # Diversity (H) of interactions 
      # diversity of interactions received 
      res[i, 1] <-  -sum(sapply(j, function(x) {
        div <- abs(imat[i, x, n]) / sum(abs(imat[i, , n]))
        div*log(div, 2)
      }))
      # diversity of interactions emitted 
      res[i, 2] <-  -sum(sapply(j, function(x) {
        div <- abs(imat[x, i, n]) / sum(abs(imat[ , i, n]))
        div*log(div, 2)
      }))
      
      # % of exploited interactions (-/+)
      vec <- sapply(1:Sp, function(j){ 
        ifelse(imat[i, j, n] > 0, 0, ifelse(  # if a_ij is competitive
          imat[j, i, n] < 0, 0, 1)) # and a_ji is facilitative
      })
      vec <- vec[-i] # remove intraspecific interaction
      res[i, 3] <- sum(vec)/length(vec) # 
      
      # % of exploiter interactions (+/-)
      vec <- sapply(1:Sp, function(j){ 
        ifelse(imat[i, j, n] < 0, 0, ifelse(  # if a_ij is facilitative
          imat[j, i, n] > 0, 0, 1)) # and a_ji is competitive
      })
      vec <- vec[-i] # remove intraspecific interaction
      res[i, 4] <- sum(vec)/length(vec) # 
      
      # % of true competitive interactions (-/-)
      vec <- sapply(1:Sp, function(j){ 
        ifelse(imat[i, j, n] > 0, # if a_ij is competitive...
           ifelse(imat[j, i, n] > 0, 1, 0), # and a_ji is competitive too
           0)
      })
      vec <- vec[-i] # remove intraspecific interaction
      res[i, 5] <- sum(vec)/length(vec)
      
      # % of co-operative interactions (+/+)
      vec <- sapply(1:Sp, function(j){ 
        ifelse(imat[i, j, n] < 0, # if a_ij is facilitative...
               ifelse(imat[j, i, n] < 0, 1, 0), # and a_ji is facilitative too
               0)
      })
      vec <- vec[-i] # remove intraspecific interaction
      res[i, 6] <- sum(vec)/length(vec)
      
    }

  return(res)
  }, simplify = 'array')
  
  metrics[ , 'H.in', ] <-  rbind(iemat[ , 1, ], matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'H.out', ] <-  rbind(iemat[ , 2, ], matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'perc.exploited', ] <-  rbind(iemat[ , 3, ], matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'perc.exploiter', ] <-  rbind(iemat[ , 4, ], matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'perc.fullcomp', ] <-  rbind(iemat[ , 5, ], matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'perc.coop', ] <-  rbind(iemat[ , 6, ], matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  
  
  # Competitive rank (how many species do you outcompete in a binary 
  # competitive outcomes matrix) 
  
  cr <- apply(imat, 3, function(cmat) { 
    # binary competitive outcomes matrix between focals 
    binmat <- sapply(1:Sp, function(i) { sapply(1:Sp, function(j) {
      ifelse(cmat[i, j] > cmat[j, i], 1, 0) # if this looks wrong it's because apply 
    })  })                                  # kinda switches around the matrix? It's weird but trust me.
    # 1 = row species outcompetes column species
    # rowsums of that (= how many species do you outcompete)
    # divide by total number of species for relative comp. rank
    rowSums(binmat)/Sp
  })
  
  metrics[ , 'comp.rank', ] <- rbind(cr, 
                                     matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  
  # Extract the response and effect parameters too 
  response <- read.csv(paste0('2.case_study/model/output/response_samples.csv'),  stringsAsFactors = F)
  effect <- read.csv(paste0('2.case_study/model/output/effect_samples.csv'), stringsAsFactors = F)
  colnames(response) <- read.csv(paste0('2.case_study/data/key_speciesID', com, '.csv'),stringsAsFactors = F)[ , 1]
  colnames(effect) <- read.csv(paste0('2.case_study/data/key_neighbourID', com, '.csv'), stringsAsFactors = F)[ , 1]
  effect <- effect[ , colnames(effect) %in% colnames(response)]
  
  metrics[ , 'response', ] <- rbind(t(response), matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'effect', ] <- t(effect)

  # sum of competition received (alphas * neighbour abundances)
  efcc <- efc(imat, com, intras)
  metrics[ , 'eff.comp', ] <- rbind(efcc$eff.comp, matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  metrics[ , 'eff.comp.intra', ] <- rbind(efcc$eff.comp.intra, matrix(data = NA, nrow = (Nb-Sp), ncol = Sm))
  
  # COMMUNITY-LEVEL ABUNDANCES
  others <- sum(abund[!abund$species %in% dimnames(metrics)[['species']], 'Num_indivs'])
  named <- abund[abund$species %in% dimnames(metrics)[['species']], c('species', 'Num_indivs')]
  total <- as.numeric(others + sum(named$Num_indivs)) # get total abundance
  named <- split(named, as.factor(named$species))
  named <- lapply(named, function(x){rep(sum(x[, 'Num_indivs']), Sm)}) # repeat by the number of samples
  named <- do.call(rbind, named)
  abunds <- rbind(named, rep(others, Sm))
  dimnames(abunds)[[1]][Nb+1] <- 'others'
  # order
  abunds <- abunds[dimnames(metrics)[['species']], ]
  # join up
  metrics[ , 'com.abund', ] <- as.numeric(abunds)
  # and get relative abundance
  metrics[ , 'rel.abund', ] <- as.numeric(abunds)/total
  
  # Get community ! 
  
  metrics[ , 'com', ] <- matrix(data = rep(as.numeric(com), Nb*Sm), nrow = Nb, ncol = Sm)
  
  # LA FIN 
  return(metrics)
}







