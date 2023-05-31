
# Organization of data in the repository:

> Simulated data are organized in different directories:

> 	* Xs: This directory contains data for predictor variables (repeated exposures data). Within this directory, the user will find an RData file under the name "resu.sim.dataX.i.RData". This file contains a list of 100 elements, each of them corresponding to a sublist for each simulated dataset. Each dataset is composed of 1200 individuals with data for 100 unique exposures and 5 time points (500 variables in total). The vector of Rho values employed for the simulation of the repeated exposures (one for each unique exposure) can also be found for each simulated dataset. These simulated datasets can be re-generated running the script (./src/scripts/dataset_simulation/XXXX.r)
	
>	* Ys: This directory contains data for outcomes. Within this directory, the user will find a sub-directory structure containing the outcome data for each of simulated scenario. Sub-directory names follow the structure "dataNUMBERexpNUMBER". In this nomenclature, data1 refers to the first scenario ("all the 5 time points of each of the true exposures are truly associated with Y"), data2 refers to the second scenario ("only a single time point of each of the true exposures is truly associated with Y"); while exp3, exp5 and exp10 refers to the number of true exposures associated with the outcome. These simulated data can be re-generated running the script (./src/scripts/dataset_simulation/YYYY.r)



