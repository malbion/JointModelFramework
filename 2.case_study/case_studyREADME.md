# Joint Model Framework / 2.case_study

This folder contains the data, model output and certain specific functions for the joint model framework applied to a case study of Western Australia annual wildflowers, as described in the associated manuscript. The code presented here is not intended to be used for other datasets, but made available in the interest of open access and reproducibility.

### How to use it: 
Use the **model.R** script - it can be run straight from the terminal, or through R (instructions in the comments at the top of the model.R file). Please note that running the model may take considerable time and processing power, and this script includes the saving of large .Rdata objects, so consider yourself forwarned! 
**model.R** will load the case study data, prepare it for the model, select identifiable interaction parameters, run the joint model, verify model convergence, save model output and selected parameters, and scale the returned interaction parameters according to a population dynamic model for annual plants with a seed bank.  
NB: **model.R** calls upon the ../1.code/**data_prep.R** function to format the data and select identifiable parameters, and the ../1.code/**joint_model.stan** file to run the model.  

### File inventory

#### data
Contains all the data files for the case study data of an annual wildflower system, called in **model.R**. Observations of seed production and neighbour abundances are stored in **fecundities0.csv**, focal species abundances are stored in **plot_species_abundances.csv,** and the demographic rates used in the population dynamic model for the scaling of interaction parameters are stored in **seed_rates.csv**. 

#### functions
Functions specific to the case study that are called in **model.R**.  
**stan_modelcheck_rem.R** contains various diagnostic functions for evaluating model fit and convergence, including the creation of traceplots which are saved in the model folder below.   
**scale_interactions.R** contains a function which scales the interaction parameters returned by the joint model framework into interaction effects which are comparable across species, as described in the Methods sections 2.6 and 2.7 of the manuscript, and in further detail in the Supplementary Methods. This scaling is specific to the system and model of plant population dynamics described in the manuscript, and should not be used for other datasets.  

#### model
Folder in which all the outputs from **model.R** and the joint model framework as applied to the case study are stored. This includes parameter samples (**model/output/**), growth rates and scaled interaction estimates (**model/transformed/**), and figures used to check and validate model convergence (**model/validation/**).  Some of these files (images, and large .Rdata files) are not version controlled on github but your own local versions can be created by running **model.R**.  



