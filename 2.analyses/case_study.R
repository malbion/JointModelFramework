

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

# Do stuff on the alphas


