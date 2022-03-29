# Joint Model Framework / 3.analyses


This folder contains the R scripts and functions used to analyse the  case study model output and create the figures for the manuscript. 

### File Inventory

#### calc_sp_ntwk_metrics.R 
Function to summarise certain metrics at a species-level from the case study results e.g. in-strength and out-strength, sum of facilitative effects, etc. These results do not appear in the manuscript. 

#### case_study.R 
Script to calculate percentage of facilitative and competitive interactions, and summarise certain interaction metrics using the calc_sp_ntwk_metrics() function.

#### figures_mss.R
Script to plot Figure 1 to 3 in the main text, and Supplementary Figures 1 and 2. 

#### ifm_vs_rem.R 
Script to plot Supplementary Figures 3, 4 and 5 aka network graphs of the NDDM (identifiable) and IFM (non-identifiable) interactions using the qgraph package. 

#### smm.csv 
Output from the calc_sp_ntwk_metrics() function, this does not appear in the manuscript.