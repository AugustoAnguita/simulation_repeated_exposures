
# Organization of software in the repository:

> All scripts employed for the generation of manuscript results are organized as follows:

> * Directory "scripts": Within this directory all employed codes are further organized into three sub-directories:
		- "dataset_simulation": This directory contains the required scripts for simulating both repeated exposure and outcome data.
		
		- "data_analysis" : This directory contains the required scripts for running all assessed techniques on the simulated datasets. Six scripts can be found within this directory, each of them prepared for running the methods in a particular scenario. Scripts names follow the structure "script_dataNUMBERexpNUMBER.R". In this nomenclature, data1 refers to the first scenario ("all the 5 time points of each of the true exposures are truly associated with Y"), data2 refers to the second scenario ("only a single time point of each of the true exposures is truly associated with Y"); while exp3, exp5 and exp10 refers to the number of true exposures associated with the outcome.
		
		- "findings_interpretation": This directory contains the required scripts for extracting the performance metrics of each method over the ...

	
> * Directory "source_functions": This directory contain a source script with the codes for running all compared methods. This script is internally called in the scripts located in the "scripts" directory, so the user will not need to modify or use it, unless is interested in introducing changes in the way methods are applied, or in case is interested in applying them to a dataset different from the simulated data.

> All required libraries are defined at the beginning of each script. We recommend the user to make sure all packages are installed in R in order to be able to run the analysis.
