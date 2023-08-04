
# Organization of results in the repository:

> Results are divided into two sub-directories: one-step and two-step approaches, each of them being divided into two sub-directories corresponding to the scenarios, each of them being divided into three sub-directories corresponding to the number of true predictors: exp3, exp5, and exp10
>	- one_step
>		* dataY1andX 
>			** exp3
>				*** 1 list for each method and reporting the result of the selection: RES1.Ewas.i.exp3.RData, RES1.DSA.i.exp3.RData, ... 
>				*** 1 list (RES1.i.all.exp3.RData) and 1 df (RES1.all.exp3.RData) combining the 6 lists mentioned above.  
>				*** the performance of each method: separetely for each simulated dataset (performance_detail_data1_exp3_XXX.Rdata) and averaged across each simulated datasets (performance_summary_data1_exp3_XXX.Rdata). When the performance are calculate to identify the true exposure at the true time point, the file name ends by "_500" while the performances calculate independently of the true time point ends by "_100".  These results can be generated using the script "./src/scripts/data_analysis/script_data1exp3.R". 
>			** exp5
>			** exp10 
		* dataY2andX (name of the resulting objects: RES2 or data2)
>			** exp3
>			** exp5
>			** exp10 
>	- two_step: the structure is similar to "one_step". The resulting object include in their name 'av' (eg., RES1av.Ewas.i.exp3.RData and performance_detail_data1av_exp3_XXX.Rdata)
>		* dataY1andX 
>			** exp3
>			** exp5
>			** exp10 
		* dataY2andX
>			** exp3
>			** exp5
>			** exp10 
