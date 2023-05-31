#######################################################################################################################################################################################
##
##
## 						Analysis of scenario 2: Only a single time point of each of the true exposures is truly associated with Y.
##                                                              **Number of true exposures associated with the outcome is 3**
##
##
## Last update: 09-05-2023
## Authors: Charline Warembourg <charline.warembourg@univ-rennes1.fr> Augusto Anguita <augusto.anguita@isglobal.org>  Xavier Basagaña <xavier.basagana@isglobal.org> 
#######################################################################################################################################################################################





########################################################
## Set up R environment 			        ##
########################################################

  
rm(list=ls())

list.of.packages <- c(
  "Matrix","MASS","parallel","glmnet","spls","MXM","DSA","dlnm","splines","mgcv","sNPLS","plyr","stringr","foreach","doParallel","ranger","palmerpenguins","tidyverse","kableExtra","future"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
for(package.i in list.of.packages){ #loading packages
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

options(max.print=100000)

print("Step 1 - Required R libraries successfully loaded")





########################################################
## Defining working directory			        ##
########################################################


setwd(here::here())
getwd()





########################################################
## Loading working datasets and source functions      ##
########################################################

source(file="./src/source_functions/functions.R")

# define number of simulated dataset (X and Y) and number of true X predictor
nsim <- 100
nexp <- 3

# load X data
load(file="./data/Xs/resu.sim.dataX.i.RData")

# load Y data
load(file="./data/Ys/dataY2andX/exp3/resu.sim.data2.exp3.RData")





########################################################
## Setting up the parallelization process             ##
########################################################

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster( # create the cluster
  n.cores, 
  type = "PSOCK"
)
print(my.cluster) # check cluster definition (optional)
doParallel::registerDoParallel(cl = my.cluster) # register it to be used by %dopar%
foreach::getDoParRegistered() # check if it is registered (optional)
foreach::getDoParWorkers() # how many workers are available? (optional)

print("Step 3 - Set up for parallelization completed")





########################################################
## Run models (ONE-STEP APPROACH)		        ##
########################################################




#####  ****************  ExWAS and ExWAS.MLR ****************  

RES2.Ewas.i.exp3<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES2.Ewas.i.exp3[[i]] = applyEwas(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
	  RES2.Ewas.i.exp3[[i]]$numsim<-i
}

save(RES2.Ewas.i.exp3,file="./results/one_step/dataY2andX/exp3/RES2.Ewas.i.exp3.RData")
load(file="./results/one_step/dataY2andX/exp3/RES2.Ewas.i.exp3.RData")




#####  ****************  ElasticNet (ENET)  **************** 

init<-RES2.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES2.Enet.i.exp3<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES2.Enet.i.exp3[[i]] = applyEnet(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
	  RES2.Enet.i.exp3[[i]]$numsim<-i
}

save(RES2.Enet.i.exp3,file="./results/one_step/dataY2andX/exp3/RES2.Enet.i.exp3.RData")




#####  ****************  sPLS  **************** 

init<-RES2.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES2.sPLS.i.exp3<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES2.sPLS.i.exp3[[i]] = applySPLS(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
	  RES2.sPLS.i.exp3[[i]]$numsim<-i
}

save(RES2.sPLS.i.exp3,file="./results/one_step/dataY2andX/exp3/RES2.sPLS.i.exp3.RData")




#####  ****************  DSA  ****************

init<-RES2.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES2.DSA.i.exp3<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES2.DSA.i.exp3[[i]] = applyDSA(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,maxsize=25)
	  RES2.DSA.i.exp3[[i]]$numsim<-i
}

save(RES2.DSA.i.exp3,file="./results/one_step/dataY2andX/exp3/RES2.DSA.i.exp3.RData")




#####  ****************  Penalized DLNM with variable selection  ****************


RES2.DLNMselect.i.exp3 <- vector("list", nsim)

RES2.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
	 
	  data.frame(applyDLNMselect(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X),i) # Apply penalized DLNM with variable selection
	  
  }
    
save(RES2.DLNMselect.i.exp3,file="./results/one_step/dataY2andX/exp3/RES2.DLNMselect.i.exp3.RData")


for(i in 1:nsim) { # this for loop allows to check the significance of the overall p-value (crossbasis) and the confidence interval of each lag : if overall p-value is significant 
	  	   # (with or without correction for multiple testing) and the CI doesn't include 0, then the exposure is considered as "selected" otherwise not.
	  RES2.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,100)
	  RES2.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES2.DLNMselect.i.exp3[[i]]$expo_name), RES2.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
	  RES2.DLNMselect.i.exp3[[i]] <- RES2.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES2.DLNMselect.i.exp3[[i]]$var))),]
	  RES2.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES2.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data2.exp3[[i]]$true.pred,1,0)
	  list.pval<-RES2.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES2.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)] # store p-values (by step of 5 = overall p-value of the crossbasis)
	  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
	  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
	  # CI includes 0 or not rep(list.pval,each=5)
	  RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES2.DLNMselect.i.exp3[[i]]$CI.inf, RES2.DLNMselect.i.exp3[[i]]$CI.sup))
	  RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
	  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
	  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.0005,1,0))
	  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
	  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
	  RES2.DLNMselect.i.exp3.v2[[i]]$DLNMpen.TP_dic <- ifelse(RES2.DLNMselect.i.exp3.v2[[i]]$p.value<0.05,1,0)
	  
}

RES2.DLNMselect.i.exp3.v2 <- RES2.DLNMselect.i.exp3

save(RES2.DLNMselect.i.exp3.v2,file="./results/one_step/dataY2andX/exp3/RES2.DLNMselect.i.exp3.v2.RData")




#####  ****************  sNPLS  ****************


future::plan(cluster, workers=my.cluster)
RES2.sNPLS.i.exp3 <- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES2.sNPLS.i.exp3[[i]] = applysNPLS(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES2.sNPLS.i.exp3[[i]]$numsim <- i
}
save(RES2.sNPLS.i.exp3,file="./results/one_step/dataY2andX/exp3/RES2.sNPLS.i.exp3.RData")





########################################################
## Run models (TWO-STEP APPROACH)		        ##
########################################################


#### MODELS STARTED WITH VARIABLES SELECTED IN ExWAS Averaged UNCORRECTED APPROACH:

#### Filtering input variables to those selected by AVG ExWAS

Selected_Xs <- vector("list", nsim)

for (i in 1:nsim) {
  
	init_varlist <- gsub("_mean","_",RES2av.Ewas.i.exp3.100[[i]]$var[which(RES2av.Ewas.i.exp3.100[[i]]$EWAS.TP.none==1)])
	
	Selected_Xs[[i]] <- resu.sim.dataX.i[[i]]$X[,unlist(lapply(as.character(init_varlist),function(x){grep(x, gsub("\\..*","_",colnames(resu.sim.dataX.i[[i]]$X)) )}))]

	}




