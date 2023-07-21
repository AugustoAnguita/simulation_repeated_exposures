#######################################################################################################################################################################################
##
##
## 								Analysis of scenario 1: All time points are truly associated with Y.
##                                                              **Number of true exposures associated with the outcome is 5**
##
##
## Last update: 21-07-2023
## Authors: Charline Warembourg <charline.warembourg@inserm.fr> Augusto Anguita <augusto.anguita@isglobal.org>  Xavier Basagaña <xavier.basagana@isglobal.org> 
#######################################################################################################################################################################################





########################################################
## Set up R environment 			      ##
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
## Defining working directory			      ##
########################################################


setwd(here::here())
getwd()




#   ********************************  ONE-STEP APPROACH   ********************************  #


########################################################
## Loading working datasets and source functions      ##
########################################################

source(file="./src/source_functions/functions.R")

# define number of simulated dataset (X and Y) and number of true X predictor
nsim <- 100
nexp <- 5

# load X data
load(file="./data/Xs/resu.sim.dataX.i.RData")

# load Y data
load(file="./data/Ys/dataY1andX/exp5/resu.sim.data1.exp5.RData")



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
## Run models (ONE-STEP APPROACH)		      ##
########################################################




#####  ****************  ExWAS and ExWAS.MLR ****************  

RES1.Ewas.i.exp5<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES1.Ewas.i.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i[[i]]$X)
	  RES1.Ewas.i.exp5[[i]]$numsim<-i
}

save(RES1.Ewas.i.exp5,file="./results/one_step/dataY1andX/exp5/RES1.Ewas.i.exp5.RData")




#####  ****************  ElasticNet (ENET)  **************** 

init<-RES1.Ewas.i.exp5[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.Enet.i.exp5<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES1.Enet.i.exp5[[i]] = applyEnet(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i[[i]]$X)
	  RES1.Enet.i.exp5[[i]]$numsim<-i
}

save(RES1.Enet.i.exp5,file="./results/one_step/dataY1andX/exp5/RES1.Enet.i.exp5.RData")




#####  ****************  sPLS  **************** 

init<-RES1.Ewas.i.exp5[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.sPLS.i.exp5<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES1.sPLS.i.exp5[[i]] = applySPLS(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i[[i]]$X)
	  RES1.sPLS.i.exp5[[i]]$numsim<-i
}

save(RES1.sPLS.i.exp5,file="./results/one_step/dataY1andX/exp5/RES1.sPLS.i.exp5.RData")




#####  ****************  DSA  ****************

init<-RES1.Ewas.i.exp5[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.DSA.i.exp5<- vector("list", nsim)

for(i in 1:nsim) {
	  print(i)
	  RES1.DSA.i.exp5[[i]] = applyDSA(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i[[i]]$X,maxsize=25)
	  RES1.DSA.i.exp5[[i]]$numsim<-i
}

save(RES1.DSA.i.exp5,file="./results/one_step/dataY1andX/exp5/RES1.DSA.i.exp5.RData")




#####  ****************  Penalized DLNM with variable selection  ****************


RES1.DLNMselect.i.exp5 <- vector("list", nsim)

RES1.DLNMselect.i.exp5 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
	 
	  data.frame(applyDLNMselect(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i[[i]]$X),i) # Apply penalized DLNM with variable selection
	  
  }
    
save(RES1.DLNMselect.i.exp5,file="./results/one_step/dataY1andX/exp5/RES1.DLNMselect.i.exp5.RData")


# The loop below allows to check the significance of the overall p-value (crossbasis) and the confidence interval of each lag : if overall p-value is significant 
# (with or without correction for multiple testing) and the CI doesn't include 0, then the exposure is considered as "selected" otherwise not.

for(i in 1:nsim) { 
	  RES1.DLNMselect.i.exp5[[i]]$num_time<-rep(1:5,100)
	  RES1.DLNMselect.i.exp5[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp5[[i]]$expo_name), RES1.DLNMselect.i.exp5[[i]]$num_time,sep = ".")
	  RES1.DLNMselect.i.exp5[[i]] <- RES1.DLNMselect.i.exp5[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp5[[i]]$var))),]
	  RES1.DLNMselect.i.exp5[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp5[[i]]$var %in% resu.sim.data1.exp5[[i]]$true.pred,1,0)
	  list.pval<-RES1.DLNMselect.i.exp5[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp5[[i]]$var,fixed = TRUE)] # store p-values (by step of 5 = overall p-value of the crossbasis)
	  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
	  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
	  # CI includes 0 or not rep(list.pval,each=5)
	  RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp5[[i]]$CI.inf, RES1.DLNMselect.i.exp5[[i]]$CI.sup))
	  RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP==FALSE,0,1)
	  RES1.DLNMselect.i.exp5[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
	  RES1.DLNMselect.i.exp5[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.0005,1,0))
	  RES1.DLNMselect.i.exp5[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
	  RES1.DLNMselect.i.exp5[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp5[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
	  RES1.DLNMselect.i.exp5.v2[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp5.v2[[i]]$p.value<0.05,1,0)
	  
}

RES1.DLNMselect.i.exp5.v2 <- RES1.DLNMselect.i.exp5

save(RES1.DLNMselect.i.exp5.v2,file="./results/one_step/dataY1andX/exp5/RES1.DLNMselect.i.exp5.v2.RData")




#####  ****************  sNPLS  ****************


future::plan(cluster, workers=my.cluster)
RES1.sNPLS.i.exp5 <- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.sNPLS.i.exp5[[i]] = applysNPLS(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES1.sNPLS.i.exp5[[i]]$numsim <- i
}
save(RES1.sNPLS.i.exp5,file="./results/one_step/dataY1andX/exp5/RES1.sNPLS.i.exp5.RData")




#####  ****************  Merge all results ****************



RES1.i.all.exp5<-vector("list",nsim)
for(i in 1:nsim) {
  RES1.i.all.exp5[[i]]<-Reduce(function(x, y) merge(x, y,by="var", all=TRUE), 
                                list(RES1.Ewas.i.exp5[[i]], RES1.Enet.i.exp5[[i]][,c("var","ENET_MIN.TP","ENET_OPT.TP")], RES1.sPLS.i.exp5[[i]][,c("var","SPLS_MIN.TP")],
                                     RES1.DSA.i.exp5[[i]][,c("var","DSA.TP")],RES1.DLNMselect.i.exp5.v2[[i]][,c("var","DLNMselect.by")],RES1.sNPLS.i.exp5[[i]][,c("var","sNPLS.TP")]))
}
save(RES1.i.all.exp5,file="RES1.i.all.exp5.RData")

RES1.all.exp5<-matrix(ncol=dim(RES1.i.all.exp5[[1]])[2])
colnames(RES1.all.exp5)<-colnames(RES1.i.all.exp5[[1]])
for (i in 1:nsim){
  RES1.all.exp5<-rbind(RES1.all.exp5,RES1.i.all.exp5[[i]])
}

RES1.all.exp5<-RES1.all.exp5[-1,c("var","numsim","true.pred","EWAS.TP.none","EWAS_LM.TP.none",
                                    "EWAS.TP.bon","EWAS_LM.TP.bon","EWAS.TP.bh","EWAS_LM.TP.bh",  
                                    "EWAS.TP.by","EWAS_LM.TP.by","ENET_MIN.TP","ENET_OPT.TP",
                                    "SPLS_MIN.TP","DSA.TP","DLNMselect.by","sNPLS.TP")]
save(RES1.all.exp5,file="RES1.all.exp5.RData")





#####  ****************  Calculate the performances (sensitivity and FDR) ****************


# Call the calculateMetrics function to identify the true exposures at the true time point (bytime=1)
results_500  <- calculateMetrics(data=RES1.all.exp5,bytime=1)

# Access the performances of each method by simulated dataset 
performance_detail_data1_exp5_500 <- results_500$performance_detail

# Access the overall performances of each method (mean and sd of all simulated datasets)
performance_summary_data1_exp5_500 <- results_500$performance_summary



# Call the calculateMetrics function to identify the true exposures independently of the true time point (bytime=0)
results_100  <- calculateMetrics(data=RES1.all.exp5,bytime=0)

# Access the performances of each method by simulated dataset 
performance_detail_data1_exp5_100 <- results_100$performance_detail

# Access the overall performances of each method (mean and sd of all simulated datasets)
performance_summary_data1_exp5_100 <- results_100$performance_summary


save(performance_detail_data1_exp5_500,file="./results/one_step/dataY1andX/exp5/performance_detail_data1_exp5_500.Rdata")
save(performance_summary_data1_exp5_500,file="./results/one_step/dataY1andX/exp5/performance_summary_data1_exp5_500.Rdata")
save(performance_detail_data1_exp5_100,file="./results/one_step/dataY1andX/exp5/performance_detail_data1_exp5_100.Rdata")
save(performance_summary_data1_exp5_100,file="./results/one_step/dataY1andX/exp5/performance_summary_data1_exp5_100.Rdata")









#   ********************************  TWO-STEP APPROACH   ********************************  #




########################################################
## Loading working datasets 			      ##
########################################################

# load X data - Averaged
load(file="./data/Xs/Xave.i.RData")

# load Y data - denominator = 100 "_mean" true predictors
load(file="./data/Ys/dataY1andX/exp5/resu.sim.data1av.exp5.100.RData")


########################################################
## Run models (TWO-STEP APPROACH)		      ##
## Step 1: selection based on averaged exposure level ##
########################################################


#####  ****************  EWAS  ****************  

RES1av.Ewas.i.exp5.100<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.Ewas.i.exp5.100[[i]] = applyEwas(data.Y=resu.sim.data1av.exp5.100[[i]],data.X=Xave.i[[i]])
  RES1av.Ewas.i.exp5.100[[i]]$numsim<-i
}

save(RES1av.Ewas.i.exp5.100,file="./results/two_step/dataY1andX/exp5/RES1av.Ewas.i.exp5.100.RData")



#####   ****************  ENET  ****************  

init<-RES1av.Ewas.i.exp5.100[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1av.Enet.i.exp5.100<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.Enet.i.exp5.100[[i]] = applyEnet(data.Y=resu.sim.data1av.exp5.100[[i]],data.X=Xave.i[[i]])
  RES1av.Enet.i.exp5.100[[i]]$numsim<-i
}

save(RES1av.Enet.i.exp5.100,file="./results/two_step/dataY1andX/exp5/RES1av.Enet.i.exp5.100.RData")


#####   ****************  sPLS  ****************  

init<-RES1av.Ewas.i.exp5.100[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1av.sPLS.i.exp5.100<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.sPLS.i.exp5.100[[i]] = applySPLS(data.Y=resu.sim.data1av.exp5.100[[i]],data.X=Xave.i[[i]])
  RES1av.sPLS.i.exp5.100[[i]]$numsim<-i
}

save(RES1av.sPLS.i.exp5.100,file="./results/two_step/dataY1andX/exp5/RES1av.sPLS.i.exp5.100.RData")


#####   ****************  DSA  ****************  

init<-RES1av.Ewas.i.exp5.100[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1av.DSA.i.exp5.100<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.DSA.i.exp5.100[[i]] = applyDSA(data.Y=resu.sim.data1av.exp5.100[[i]],data.X=Xave.i[[i]],maxsize = 25)
  RES1av.DSA.i.exp5.100[[i]]$numsim<-i
}

save(RES1av.DSA.i.exp5.100,file="./results/two_step/dataY1andX/exp5/RES1av.DSA.i.exp5.100.RData")



########################################################
## Run models (TWO-STEP APPROACH)		      ##
## Step 2: Selection based on time-specific exposure  ## 
##         level following selection at step 1 	      ##
########################################################

#### Reduced number of X variables: 
# 1) identify the averaged data significant at the first step
# 2) reduced X data (time-specific variables) to significant X averaged 
# 3) Apply again the method to these reduced set of X data  (time-specific variables) 



#####   ****************  ExWAS  ****************

RES1av.Ewas.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.Ewas.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.Ewas.i.exp5.v2[[i]]$var<-substr((RES1av.Ewas.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.Ewas.i.exp5.v2[[i]]$var))-5)
  RES1av.Ewas.i.exp5.v3[[i]]<-RES1av.Ewas.i.exp5.v2[[i]][RES1av.Ewas.i.exp5.v2[[i]]$EWAS.TP.none==1,]
  var.sign.exp5[[i]]<-RES1av.Ewas.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.Ewas.i.redExWAS.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.Ewas.i.redExWAS.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.Ewas.i.redExWAS.exp5[[i]]$numsim<-i
  RES1av.Ewas.i.redExWAS.exp5[[i]]$candidate.red[RES1av.Ewas.i.redExWAS.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.Ewas.i.redExWAS.exp5[[i]]<-RES1av.Ewas.i.redExWAS.exp5[[i]][colnames(RES1av.Ewas.i.redExWAS.exp5[[i]])%in%c("var","val","Est","pVal","EWAS.none","EWAS.TP.none","true.pred","nsim","candidate.red")]
}

save(RES1av.Ewas.i.redExWAS.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.Ewas.i.redExWAS.exp5.RData")

#####   ****************  ExWAS - LM  ****************

RES1av.EwasLM.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasLM.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasLM.i.exp5.v2[[i]]$var<-substr((RES1av.EwasLM.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasLM.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasLM.i.exp5.v3[[i]]<-RES1av.EwasLM.i.exp5.v2[[i]][RES1av.EwasLM.i.exp5.v2[[i]]$EWAS_LM.TP.none==1,]
  var.sign.exp5[[i]]<-RES1av.EwasLM.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasLM.i.redExWASLM.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasLM.i.redExWASLM.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasLM.i.redExWASLM.exp5[[i]]$numsim<-i
  RES1av.EwasLM.i.redExWASLM.exp5[[i]]$candidate.red[RES1av.EwasLM.i.redExWASLM.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasLM.i.redExWASLM.exp5[[i]]<-RES1av.EwasLM.i.redExWASLM.exp5[[i]][colnames(RES1av.EwasLM.i.redExWASLM.exp5[[i]])%in%c("var","val","Est","pVal","EWAS_LM.none","EWAS_LM.TP.none","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasLM.i.redExWASLM.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasLM.i.redExWASLM.exp5.RData")

#####   ****************  ExWAS - Bonferroni  ****************

RES1av.EwasBonf.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasBonf.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasBonf.i.exp5.v2[[i]]$var<-substr((RES1av.EwasBonf.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasBonf.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasBonf.i.exp5.v3[[i]]<-RES1av.EwasBonf.i.exp5.v2[[i]][RES1av.EwasBonf.i.exp5.v2[[i]]$EWAS.TP.bon==1,]
  var.sign.exp5[[i]]<-RES1av.EwasBonf.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasBonf.i.redExWASBonf.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasBonf.i.redExWASBonf.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasBonf.i.redExWASBonf.exp5[[i]]$numsim<-i
  RES1av.EwasBonf.i.redExWASBonf.exp5[[i]]$candidate.red[RES1av.EwasBonf.i.redExWASBonf.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasBonf.i.redExWASBonf.exp5[[i]]<-RES1av.EwasBonf.i.redExWASBonf.exp5[[i]][colnames(RES1av.EwasBonf.i.redExWASBonf.exp5[[i]])%in%c("var","val","Est","pVal","EWAS.TP.bon","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasBonf.i.redExWASBonf.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasBonf.i.redExWASBonf.exp5.RData")


#####   ****************  ExWAS - Bonf LM  ****************

RES1av.EwasBonfLM.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasBonfLM.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasBonfLM.i.exp5.v2[[i]]$var<-substr((RES1av.EwasBonfLM.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasBonfLM.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasBonfLM.i.exp5.v3[[i]]<-RES1av.EwasBonfLM.i.exp5.v2[[i]][RES1av.EwasBonfLM.i.exp5.v2[[i]]$EWAS_LM.TP.bon==1,]
  var.sign.exp5[[i]]<-RES1av.EwasBonfLM.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasBonfLM.i.redExWASBonfLM.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]]$numsim<-i
  RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]]$candidate.red[RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]]<-RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]][colnames(RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]])%in%c("var","val","Est","pVal","EWAS_LM.bon","EWAS_LM.TP.bon","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasBonfLM.i.redExWASBonfLM.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasBonfLM.i.redExWASBonfLM.exp5.RData")



#####   ****************  ExWAS - BH  ****************

RES1av.EwasBH.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasBH.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasBH.i.exp5.v2[[i]]$var<-substr((RES1av.EwasBH.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasBH.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasBH.i.exp5.v3[[i]]<-RES1av.EwasBH.i.exp5.v2[[i]][RES1av.EwasBH.i.exp5.v2[[i]]$EWAS.TP.bh==1,]
  var.sign.exp5[[i]]<-RES1av.EwasBH.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasBH.i.redExWASBH.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasBH.i.redExWASBH.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasBH.i.redExWASBH.exp5[[i]]$numsim<-i
  RES1av.EwasBH.i.redExWASBH.exp5[[i]]$candidate.red[RES1av.EwasBH.i.redExWASBH.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasBH.i.redExWASBH.exp5[[i]]<-RES1av.EwasBH.i.redExWASBH.exp5[[i]][colnames(RES1av.EwasBH.i.redExWASBH.exp5[[i]])%in%c("var","val","Est","pVal","EWAS.TP.bh","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasBH.i.redExWASBH.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasBH.i.redExWASBH.exp5.RData")


#####   ****************  ExWAS - BH LM  ****************

RES1av.EwasBHLM.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasBHLM.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasBHLM.i.exp5.v2[[i]]$var<-substr((RES1av.EwasBHLM.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasBHLM.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasBHLM.i.exp5.v3[[i]]<-RES1av.EwasBHLM.i.exp5.v2[[i]][RES1av.EwasBHLM.i.exp5.v2[[i]]$EWAS_LM.TP.bh==1,]
  var.sign.exp5[[i]]<-RES1av.EwasBHLM.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasBHLM.i.redExWASBHLM.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]]$numsim<-i
  RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]]$candidate.red[RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]]<-RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]][colnames(RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]])%in%c("var","val","Est","pVal","EWAS_LM.bh","EWAS_LM.TP.bh","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasBHLM.i.redExWASBHLM.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasBHLM.i.redExWASBHLM.exp5.RData")


#####   ****************  ExWAS - BY  ****************

RES1av.EwasBY.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasBY.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasBY.i.exp5.v2[[i]]$var<-substr((RES1av.EwasBY.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasBY.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasBY.i.exp5.v3[[i]]<-RES1av.EwasBY.i.exp5.v2[[i]][RES1av.EwasBY.i.exp5.v2[[i]]$EWAS.TP.by==1,]
  var.sign.exp5[[i]]<-RES1av.EwasBY.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasBY.i.redExWASBY.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasBY.i.redExWASBY.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasBY.i.redExWASBY.exp5[[i]]$numsim<-i
  RES1av.EwasBY.i.redExWASBY.exp5[[i]]$candidate.red[RES1av.EwasBY.i.redExWASBY.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasBY.i.redExWASBY.exp5[[i]]<-RES1av.EwasBY.i.redExWASBY.exp5[[i]][colnames(RES1av.EwasBY.i.redExWASBY.exp5[[i]])%in%c("var","val","Est","pVal","EWAS.TP.by","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasBY.i.redExWASBY.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasBY.i.redExWASBY.exp5.RData")


#####   ****************  ExWAS - BY LM  ****************

RES1av.EwasBYLM.i.exp5.v2<-RES1av.Ewas.i.exp5.100
RES1av.EwasBYLM.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.EwasBYLM.i.exp5.v2[[i]]$var<-substr((RES1av.EwasBYLM.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.EwasBYLM.i.exp5.v2[[i]]$var))-5)
  RES1av.EwasBYLM.i.exp5.v3[[i]]<-RES1av.EwasBYLM.i.exp5.v2[[i]][RES1av.EwasBYLM.i.exp5.v2[[i]]$EWAS_LM.TP.by==1,]
  var.sign.exp5[[i]]<-RES1av.EwasBYLM.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.EwasBYLM.i.redExWASBYLM.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]] = applyEwas(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]]$numsim<-i
  RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]]$candidate.red[RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]]<-RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]][colnames(RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]])%in%c("var","val","Est","pVal","EWAS_LM.by","EWAS_LM.TP.by","true.pred","nsim","candidate.red")]
}

save(RES1av.EwasBYLM.i.redExWASBYLM.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.EwasBYLM.i.redExWASBYLM.exp5.RData")


#####   ****************  Enet - MIN  ****************

RES1av.Enet.i.exp5.v2<-RES1av.Enet.i.exp5.100
RES1av.Enet.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.Enet.i.exp5.v2[[i]]$var<-substr((RES1av.Enet.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.Enet.i.exp5.v2[[i]]$var))-5)
  RES1av.Enet.i.exp5.v3[[i]]<-RES1av.Enet.i.exp5.v2[[i]][RES1av.Enet.i.exp5.v2[[i]]$ENET_MIN.TP==1,]
  var.sign.exp5[[i]]<-RES1av.Enet.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.Enet.i.redEnetMin.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.Enet.i.redEnetMin.exp5[[i]] = applyEnet(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.Enet.i.redEnetMin.exp5[[i]]$numsim<-i
  RES1av.Enet.i.redEnetMin.exp5[[i]]$candidate.red[RES1av.Enet.i.redEnetMin.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
}
save(RES1av.Enet.i.redEnetMin.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.Enet.i.redEnetMin.exp5.RData")


#####   ****************  Enet - OPT  ****************

RES1av.Enet.i.exp5.v2<-RES1av.Enet.i.exp5.100
RES1av.Enet.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.Enet.i.exp5.v2[[i]]$var<-substr((RES1av.Enet.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.Enet.i.exp5.v2[[i]]$var))-5)
  RES1av.Enet.i.exp5.v3[[i]]<-RES1av.Enet.i.exp5.v2[[i]][RES1av.Enet.i.exp5.v2[[i]]$ENET_OPT.TP==1,]
  var.sign.exp5[[i]]<-RES1av.Enet.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}


RES1av.Enet.i.redEnetOpt.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  if (ncol(resu.sim.dataX.i.red[[i]]$X)==0){
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$var<-init$var
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$numsim<-rep(i,500)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$ENET_MIN<-rep(NA,500)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$ENET_OPT<-rep(NA,500)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$ENET_MIN.TP<-rep(NA,500)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$ENET_OPT.TP<-rep(0,500)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$candidate.red<-rep(0,500)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]<-as.data.frame(RES1av.Enet.i.redEnetOpt.exp5[[i]])
  }else{
    RES1av.Enet.i.redEnetOpt.exp5[[i]] = applyEnet(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$numsim<-i
    RES1av.Enet.i.redEnetOpt.exp5[[i]]$candidate.red[RES1av.Enet.i.redEnetOpt.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  }
}
save(RES1av.Enet.i.redEnetOpt.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.Enet.i.redEnetOpt.exp5.RData")



#####   ****************  sPLS   ****************

RES1av.sPLS.i.exp5.v2<-RES1av.sPLS.i.exp5.100
RES1av.sPLS.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.sPLS.i.exp5.v2[[i]]$var<-substr((RES1av.sPLS.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.sPLS.i.exp5.v2[[i]]$var))-5)
  RES1av.sPLS.i.exp5.v3[[i]]<-RES1av.sPLS.i.exp5.v2[[i]][RES1av.sPLS.i.exp5.v2[[i]]$SPLS_MIN.TP==1,]
  var.sign.exp5[[i]]<-RES1av.sPLS.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.sPLS.i.redsPLS.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1av.sPLS.i.redsPLS.exp5[[i]] = applySPLS(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X)
  RES1av.sPLS.i.redsPLS.exp5[[i]]$numsim<-i
  RES1av.sPLS.i.redsPLS.exp5[[i]]$candidate.red[RES1av.sPLS.i.redsPLS.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
}
save(RES1av.sPLS.i.redsPLS.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.sPLS.i.redsPLS.exp5.RData")


#####   ****************  DSA   ****************

RES1av.DSA.i.exp5.v2<-RES1av.DSA.i.exp5.100
RES1av.DSA.i.exp5.v3<-list("vector",nsim)
resu.sim.dataX.i.red<-resu.sim.dataX.i
var.sign.exp5<-list("vector",nsim)
for (i in 1:100){
  RES1av.DSA.i.exp5.v2[[i]]$var<-substr((RES1av.DSA.i.exp5.v2[[i]]$var),1,nchar(as.character(RES1av.DSA.i.exp5.v2[[i]]$var))-5)
  RES1av.DSA.i.exp5.v3[[i]]<-RES1av.DSA.i.exp5.v2[[i]][RES1av.DSA.i.exp5.v2[[i]]$DSA.TP==1,]
  var.sign.exp5[[i]]<-RES1av.DSA.i.exp5.v3[[i]]$var
  resu.sim.dataX.i.red[[i]]$X<-resu.sim.dataX.i.red[[i]]$X[,substr(colnames(resu.sim.dataX.i.red[[i]]$X),1,nchar(colnames(resu.sim.dataX.i.red[[i]]$X))-2)%in%(var.sign.exp5[[i]])]
}

RES1av.DSA.i.redDSA.exp5<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  if (dim(resu.sim.dataX.i.red[[i]]$X)[2]==0){
    RES1av.DSA.i.redDSA.exp5[[i]]<-NULL
  }else {
    RES1av.DSA.i.redDSA.exp5[[i]] = applyDSA(data.Y=resu.sim.data1.exp5[[i]],data.X=resu.sim.dataX.i.red[[i]]$X,maxsize = 25)
    RES1av.DSA.i.redDSA.exp5[[i]]$numsim<-i
    RES1av.DSA.i.redDSA.exp5[[i]]$candidate.red[RES1av.DSA.i.redDSA.exp5[[i]]$var %in% colnames(resu.sim.dataX.i.red[[i]]$X)]<-1 
  }
}
save(RES1av.DSA.i.redDSA.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.DSA.i.redDSA.exp5.RData")








#####  ****************  Merge all results ****************


RES1av.i.all.exp5<-vector("list",nsim)
for(i in 1:nsim) {
  RES1av.i.all.exp5[[i]]<-Reduce(function(x, y) merge(x, y,all=TRUE,by="var"), 
                      list(RES1.Ewas.i.exp5[[i]][,c("var","val","Est","pVal","true.pred","numsim")],RES1av.Ewas.i.redExWAS.exp5[[i]][,c("var","EWAS.TP.none")], RES1av.EwasLM.i.redExWASLM.exp5[[i]][,c("var","EWAS_LM.TP.none")], RES1av.EwasBonf.i.redExWASBonf.exp5[[i]][,c("var","EWAS.TP.bon")],
                           RES1av.EwasBonfLM.i.redExWASBonfLM.exp5[[i]][,c("var","EWAS_LM.TP.bon")],RES1av.EwasBH.i.redExWASBH.exp5[[i]][,c("var","EWAS.TP.bh")],RES1av.EwasBHLM.i.redExWASBHLM.exp5[[i]][,c("var","EWAS_LM.TP.bh")],
                           RES1av.EwasBY.i.redExWASBY.exp5[[i]][,c("var","EWAS.TP.by")],RES1av.EwasBYLM.i.redExWASBYLM.exp5[[i]][,c("var","EWAS_LM.TP.by")],RES1av.Enet.i.redEnetMin.exp5[[i]][,c("var","ENET_MIN.TP")],
                           RES1av.Enet.i.redEnetOpt.exp5[[i]][,c("var","ENET_OPT.TP")],RES1av.sPLS.i.redsPLS.exp5[[i]][,c("var","SPLS_MIN.TP")],RES1av.DSA.i.redDSA.exp5[[i]][,c("var","DSA.TP")]))
}


RES1av.all.exp5<-matrix(ncol=dim(RES1av.i.all.exp5[[1]])[2])
colnames(RES1av.all.exp5)<-colnames(RES1av.i.all.exp5[[1]])
for (i in 1:nsim){
  RES1av.all.exp5<-rbind(RES1av.all.exp5,RES1av.i.all.exp5[[i]])
}


RES1av.all.exp5<-RES1av.all.exp5[!is.na(RES1av.all.exp5$var),]

RES1av.all.exp5$EWAS.TP.none[is.na(RES1av.all.exp5$EWAS.TP.none)]<-"0"
RES1av.all.exp5$EWAS_LM.TP.none[is.na(RES1av.all.exp5$EWAS_LM.TP.none)]<-"0"
RES1av.all.exp5$EWAS.TP.bon[is.na(RES1av.all.exp5$EWAS.TP.bon)]<-"0"
RES1av.all.exp5$EWAS_LM.TP.bon[is.na(RES1av.all.exp5$EWAS_LM.TP.bon)]<-"0"
RES1av.all.exp5$EWAS.TP.bh[is.na(RES1av.all.exp5$EWAS.TP.bh)]<-"0"
RES1av.all.exp5$EWAS_LM.TP.bh[is.na(RES1av.all.exp5$EWAS_LM.TP.bh)]<-"0"
RES1av.all.exp5$EWAS.TP.by[is.na(RES1av.all.exp5$EWAS.TP.by)]<-"0"
RES1av.all.exp5$EWAS_LM.TP.by[is.na(RES1av.all.exp5$EWAS_LM.TP.by)]<-"0"
RES1av.all.exp5$ENET_MIN.TP[is.na(RES1av.all.exp5$ENET_MIN.TP)]<-"0"
RES1av.all.exp5$ENET_OPT.TP[is.na(RES1av.all.exp5$ENET_OPT.TP)]<-"0"
RES1av.all.exp5$SPLS_MIN.TP[is.na(RES1av.all.exp5$SPLS_MIN.TP)]<-"0"
RES1av.all.exp5$DSA.TP[is.na(RES1av.all.exp5$DSA.TP)]<-"0"
RES1av.all.exp5$true.pred[is.na(RES1av.all.exp5$true.pred)]<-"0"


save(RES1av.all.exp5,file="./results/two_step/dataY1andX/exp5/RES1av.all.exp5.RData")




#####  ****************  Calculate the performances (sensitivity and FDR) ****************


# Call the calculateMetrics function to identify the true exposures at the true time point (bytime=1)
results_500  <- calculateMetrics(data=RES1av.all.exp5,bytime=1)

# Access the performances of each method by simulated dataset 
performance_detail_data1av_exp5_500 <- results_500$performance_detail

# Access the overall performances of each method (mean and sd of all simulated datasets)
performance_summary_data1av_exp5_500 <- results_500$performance_summary



# Call the calculateMetrics function to identify the true exposures independently of the true time point (bytime=0)
results_100  <- calculateMetrics(data=RES1av.all.exp5,bytime=0)

# Access the performances of each method by simulated dataset 
performance_detail_data1av_exp5_100 <- results_100$performance_detail

# Access the overall performances of each method (mean and sd of all simulated datasets)
performance_summary_data1av_exp5_100 <- results_100$performance_summary

save(performance_detail_data1av_exp5_500,file="./results/two_step/dataY1andX/exp5/performance_detail_data1av_exp5_500.Rdata")
save(performance_summary_data1av_exp5_500,file="./results/two_step/dataY1andX/exp5/performance_summary_data1av_exp5_500.Rdata")
save(performance_detail_data1av_exp5_100,file="./results/two_step/dataY1andX/exp5/performance_detail_data1av_exp5_100.Rdata")
save(performance_summary_data1av_exp5_100,file="./results/two_step/dataY1andX/exp5/performance_summary_data1av_exp5_100.Rdata")














# /!\ Augusto's code? To kkep or not ??


#### MODELS STARTED WITH VARIABLES SELECTED IN ExWAS Averaged UNCORRECTED APPROACH:

#### Filtering input variables to those selected by AVG ExWAS

# Selected_Xs <- vector("list", nsim)

# for (i in 1:nsim) {
  
# 	init_varlist <- gsub("_mean","_",RES1av.Ewas.i.exp5.100[[i]]$var[which(RES1av.Ewas.i.exp5.100[[i]]$EWAS.TP.none==1)])
	
#	Selected_Xs[[i]] <- resu.sim.dataX.i[[i]]$X[,unlist(lapply(as.character(init_varlist),function(x){grep(x, gsub("\\..*","_",colnames(resu.sim.dataX.i[[i]]$X)) )}))]

#	}




