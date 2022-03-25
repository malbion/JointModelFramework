# Joint Model Framework / 1.code

This folder contains all the code needed to run and implement the joint model framework on a simulated data set. We recommend you use the code presented here if you would like to run the joint model framework on your own dataset. 

### How to use it: 
1- open the **master.R** script in R 
2- the script will call up the **simul_data()** function to simulate data, and the **data_prep()** function to check parameter identifiability and prepare the data for the format preferred by STAN
3- the script will then run the model using the **joint_model.stan** file
4- Tada! 


### File inventory: 
#### data_prep.R
This contains the data_prep() function which is called in master.R
data_prep() takes a dataframe object where rows are observations of individual performance and columns include the densities of potential interaction partners or neighbours. It returns the data in a list object which is the correct format with which to run the model on. This list includes Q, the matrix of identifiable interaction parameters, how this matrix is calculated can be found in this file under the heading 'Matrix of inferrable interactions'. 

#### joint_model.stan
This file describes the joint model in the STAN programing language. Comments include the corresponding names of the various parameters as named in the associated manuscript. This file can be modified as necessary to model the response variable with a different distribution and link function than the one used in our case study (negative binomial).

#### master.R
This is the 'master' script which calls all the other functions in this folder in order to simulate some data, prepare it for the model, select identifiable parameters, and then run the model in STAN. 

#### nddm.stan 
STAN model file for the NDDM only, exists for troubleshooting purposes.

#### rim.stan
STAN model file for the RIM only, exists for troubleshooting purposes.

#### simul_data.R 
This file contains the simul_data() function which simulates a simple dataset on which the joint model framework can be run.

