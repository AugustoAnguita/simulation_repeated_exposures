rm(list=ls())
gc()




################################# Set R environment 


list.of.packages <- c(
  "foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra"
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

print("Step 1 - Required R libraries successfully loaded")




################################# Set working directories

#Run R in /PROJECTES/SHARED/AugustoAnguita/ATHLETE_sim_long/DLNM/ANALYSIS/simulation_mars_2022
here::here() # check current wd // https://cran.r-project.org/web/packages/here/vignettes/here.html
here::i_am("dataY1andX/exp3/data1exp3_FINAL_MODELS_09_05_22.R") # Declare the location of the current script and set wd / The project root is initialized (one step up), consistent with the location of the current script or report
setwd(here::here()) # working directory now: "/PROJECTES/SHARED/AugustoAnguita/ATHLETE_sim_long/DLNM/ANALYSIS/simulation_mars_2022"
getwd()




################################# Load input files

source(file="Baseline/Functions method simulation.R")

# Define number of simulated datasets (X and Y) and number of true X predictors
nsim<-100
nexp<-3

# Load X data
load(file="Baseline/resu.sim.dataX.i.RData")

# Load Y data
load(file="dataY1andX/exp3/resu.sim.data1.exp3.RData")

# Load ExWAS selected variables for initialization
load(file="AVG_ExWAS_results/RES1av.Ewas.i.exp3.100.RData")

print("Step 2 - Loaded simulated predictor and outcome datasets")




################################# RUN MODELS

#### Setting up the parallelization process

parallel::detectCores()
n.cores <- 24 #parallel::detectCores() - 1
my.cluster <- parallel::makeCluster( # create the cluster
  n.cores, 
  type = "PSOCK"
)
print(my.cluster) # check cluster definition (optional)
doParallel::registerDoParallel(cl = my.cluster) # register it to be used by %dopar%
foreach::getDoParRegistered() # check if it is registered (optional)
foreach::getDoParWorkers() # how many workers are available? (optional)

print("Step 3 - Set up for parallelization completed")


######## ######## ######## ######## ######## ######## ######## ######## ######## ########
#### MODELS STARTED WITH VARIABLES SELECTED IN ExWAS Averaged UNCORRECTED APPROACH: ####
######## ######## ######## ######## ######## ######## ######## ######## ######## ########


### Filtering input variables to those selected by AVG ExWAS

Selected_Xs <- vector("list", nsim)

for (i in 1:nsim) {
  
	init_varlist <- gsub("_mean","_",RES1av.Ewas.i.exp3.100[[i]]$var[which(RES1av.Ewas.i.exp3.100[[i]]$EWAS.TP.none==1)])
	
	Selected_Xs[[i]] <- resu.sim.dataX.i[[i]]$X[,unlist(lapply(as.character(init_varlist),function(x){grep(x, gsub("\\..*","_",colnames(resu.sim.dataX.i[[i]]$X)) )}))]

	}


###########################################################
###### MODEL 1: Penalized DLM from Uncorrected ExWAS ######
###########################################################

# Check the significance of the overall p-value (crossbasis) and the confidence interval of each lag : if overall p-value is significant 
# (with or without correction for multiple testing) and the CI doesn't include 0, then considered as "selected" otherwise not (denominator=500).

RES1.DLNMpen.i.exp3_AVGEXWAS <- vector("list", nsim)

for(i in 1:nsim) {
  
  print(i)
  
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]] = applyDLNMpenAVG(data.Y=resu.sim.data1.exp3[[i]],data.X=Selected_Xs[[i]])
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$numsim <- i
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$true.pred<-ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  # P-value of the crossbasis (overall p-value)
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP<-as.numeric(as.character(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP))
  # Confidence interval of each lag includes 0 or not
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP<-as.numeric(as.character(ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP=="TRUE",1,0)))
  # Store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP[grep(".1",RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$var,fixed = TRUE)]
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.none<-ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP %in% c(0),0,ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP<0.05,1,0))
  
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.bonf<-ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP %in% c(0),0,ifelse(rep(p.adjust(list.pval,"bonferroni"),each=5)<0.05,1,0))
                                                
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.bh<- ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.by<-ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]<-RES1.DLNMpen.i.exp3_AVGEXWAS[[i]][,c("var","numsim","true.pred","DLNMpen.TP","estimate","CI.inf","CI.sup","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by","DLNMpen_lag.TP")]
  # RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$sign.lag<-ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP>0.05,0,
  #                                      ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP<0.05 & RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen_lag.TP==0,0,1))
}

#save(RES1.DLNMpen.i.exp3_AVGEXWAS,file="dataY1andX/exp3/RES1.DLNMpen.i.exp3_AVGEXWAS.RData")


#*#
#load("dataY1andX/exp3/RES1.DLNMpen.i.exp3_AVGEXWAS.RData")

for(i in 1:nsim) {
  RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMpen.i.exp3_AVGEXWAS[[i]]$DLNMpen.TP<0.05,1,0)
}

save(RES1.DLNMpen.i.exp3_AVGEXWAS,file="dataY1andX/exp3/RES1.DLNMpen.i.exp3_AVGEXWAS_FINAL.RData")
#*#

print("Step 4 - Completed Univariate model penalized DLNM without variable selection starting with Uncorrected Averaged ExWAS selected variables (output: RES1.DLNMpen.i.exp3_AVGEXWAS_FINAL.RData)")


#####################################################################################
###### MODEL 2: Penalized DLNM  with variable selection from Uncorrected ExWAS ######
#####################################################################################

RES1.DLNMselect.i.exp3.AVGEXWAS <- vector("list", nsim)

RES1.DLNMselect.i.exp3.AVGEXWAS <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
 
  data.frame(applyDLNMselectAVG(data.Y=resu.sim.data1.exp3[[i]],data.X=Selected_Xs[[i]]),i)
  
  }
    
#save(RES1.DLNMselect.i.exp3.AVGEXWAS,file="dataY1andX/exp3/RES1.DLNMselect.i.exp3.AVGEXWAS.RData")

print("Step 5 - Completed Multivariate model penalized DLNM with variable selection starting with Uncorrected Averaged ExWAS selected variables (output: RES1.DLNMselect.i.exp3.AVGEXWAS.RData)")

for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$num_time<-rep(1:5,nrow(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]])/5)
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$expo_name), RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$num_time,sep = ".")
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]] <- RES1.DLNMselect.i.exp3.AVGEXWAS[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$var))),]
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$var,fixed = TRUE)]
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)

  # CI includes 0 or not rep(list.pval,each=5)
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$CI.inf, RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
  
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.0005,1,0))
  
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  #RES1.DLNMselect.i.exp3.AVGEXWAS[[i]]<-RES1.DLNMselect.i.exp3.AVGEXWAS[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.AVGEXWAS.v2<-RES1.DLNMselect.i.exp3.AVGEXWAS

#save(RES1.DLNMselect.i.exp3.AVGEXWAS.v2,file="dataY1andX/exp3/RES1.DLNMselect.i.exp3.AVGEXWAS.v2.RData")


#*#
#load("dataY1andX/exp3/RES1.DLNMselect.i.exp3.AVGEXWAS.v2.RData")

for(i in 1:nsim) {
  RES1.DLNMselect.i.exp3.AVGEXWAS.v2[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp3.AVGEXWAS.v2[[i]]$p.value<0.05,1,0)
}

save(RES1.DLNMselect.i.exp3.AVGEXWAS.v2,file="dataY1andX/exp3/RES1.DLNMselect.i.exp3.AVGEXWAS.v2_FINAL.RData")
#*#

print("Step 6 - Completed Multivariate model penalized DLNM with variable selection with additional columns  starting with Uncorrected Averaged ExWAS selected variables  (output: RES1.DLNMselect.i.exp3.AVGEXWAS.v2_FINAL.RData)")




###########################################################################################################
###### MODEL 3: Penalized DLNM with variable selection and backward prunning from Uncorrected ExWAS ######
###########################################################################################################

padj <- "none" #This argument allows to make more difficult a predictor to stay in the model during the backward prunning (correcting pvalues by the number of predictors in the model)

RES1.DLNMselect.i.exp3 <- vector("list", nsim)

RES1.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
                                                         "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
                                                           
                                                           data.frame(applyDLNMbackwardselAVG(data.Y=resu.sim.data1.exp3[[i]],data.X=Selected_Xs[[i]],padj="none"),i)[,-10]
                                                           
                             }

#save(RES1.DLNMselect.i.exp3,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.AVGEXWAS.",as.character(padj),".RData"))

print("Step 7 - Completed model penalized DLNM with variable selection Backward prunning starting with Uncorrected Averaged ExWAS selected variables (output: RES1.DLNMselect.i.exp3.backward.AVGEXWAS...)")

for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES1.DLNMselect.i.exp3[[i]])/5)
  RES1.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3[[i]]$expo_name), RES1.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES1.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES1.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp3[[i]]$var))),]
  RES1.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3[[i]]$CI.inf, RES1.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.bonferroni<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  #  RES1.DLNMselect.i.exp3[[i]]<-RES1.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.v2_backward <- RES1.DLNMselect.i.exp3

#save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.AVGEXWAS.v2.",as.character(padj),".RData"))

#*#
#load(paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.AVGEXWAS.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES1.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp3.v2_backward[[i]]$p.value<0.05,1,0)
  
}

save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.AVGEXWAS.v2.",as.character(padj),"_FINAL.RData"))
#*#

print("Step 8 - Completed model penalized DLNM with variable selection Backward prunning with additional columns starting with Uncorrected Averaged ExWAS selected variables  (output: RES1.DLNMselect.i.exp3.backward..AVGEXWAS.v2...)")




######## ######## ######## ######## ######## ######## ######## ######## ######## ########
#### BACKWARD MODEL WITH BONF, BH, BY and none correction (to compare with DLNMsel:  ####
######## ######## ######## ######## ######## ######## ######## ######## ######## ########




###########################################################
###### MODEL 4: BACKWARD NONE CORRECTION	       ######
###########################################################


####  Penalized DLNM with variable selection and backward prunning (Xavier's proposal)

padj <- "none" #This argument allows to make more difficult a predictor to stay in the model during the backward prunning (correcting pvalues by the number of predictors in the model)

RES1.DLNMselect.i.exp3 <- vector("list", nsim)

RES1.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
  
  data.frame(applyDLNMbackwardsel(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,padj="none"),i)[,-10]
  
}

#save(RES1.DLNMselect.i.exp3,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.",as.character(padj),".RData"))

print("Step 9 - Completed model penalized DLNM with variable selection Backward prunning (NONE correction) (output: RES1.DLNMselect.i.exp3.backward...)")


for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES1.DLNMselect.i.exp3[[i]])/5)
  RES1.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3[[i]]$expo_name), RES1.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES1.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES1.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp3[[i]]$var))),]
  RES1.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3[[i]]$CI.inf, RES1.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.bonferroni<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  #  RES1.DLNMselect.i.exp3[[i]]<-RES1.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.v2_backward<-RES1.DLNMselect.i.exp3

#save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

#*#
#load(paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES1.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp3.v2_backward[[i]]$p.value<0.05,1,0)
  
}

save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),"_FINAL.RData"))
#*#

print("Step 10 - Completed model penalized DLNM with variable selection Backward prunning with additional columns (NONE correction) (output: RES1.DLNMselect.i.exp3.backward.v2...)")




###########################################################
###### MODEL 5: BACKWARD BONFERRONI CORRECTION       ######
###########################################################


####  Penalized DLNM with variable selection and backward prunning (Xavier's proposal)

padj <- "bonf" #This argument allows to make more difficult a predictor to stay in the model during the backward prunning (correcting pvalues by the number of predictors in the model)

RES1.DLNMselect.i.exp3 <- vector("list", nsim)

RES1.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
  
  data.frame(applyDLNMbackwardsel(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,padj="bonf"),i)[,-10]
  
}

#save(RES1.DLNMselect.i.exp3,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.",as.character(padj),".RData"))

print("Step 11 - Completed model penalized DLNM with variable selection Backward prunning (Bonferroni correction) (output: RES1.DLNMselect.i.exp3.backward...)")


for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES1.DLNMselect.i.exp3[[i]])/5)
  RES1.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3[[i]]$expo_name), RES1.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES1.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES1.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp3[[i]]$var))),]
  RES1.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3[[i]]$CI.inf, RES1.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.bonferroni<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  #  RES1.DLNMselect.i.exp3[[i]]<-RES1.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.v2_backward<-RES1.DLNMselect.i.exp3

#save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

#*#
#load(paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES1.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp3.v2_backward[[i]]$p.value<0.05,1,0)
  
}

save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),"_FINAL.RData"))
#*#

print("Step 12 - Completed model penalized DLNM with variable selection Backward prunning with additional columns (Bonferroni correction) (output: RES1.DLNMselect.i.exp3.backward.v2...)")




###########################################################
###### MODEL 6: BACKWARD BH CORRECTION               ######
###########################################################


####  Penalized DLNM with variable selection and backward prunning (Xavier's proposal)

padj <- "bh" #This argument allows to make more difficult a predictor to stay in the model during the backward prunning (correcting pvalues by the number of predictors in the model)

RES1.DLNMselect.i.exp3 <- vector("list", nsim)

RES1.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
  
  data.frame(applyDLNMbackwardsel(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,padj="bh"),i)[,-10]
  
}

#save(RES1.DLNMselect.i.exp3,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.",as.character(padj),".RData"))

print("Step 13 - Completed model penalized DLNM with variable selection Backward prunning (BH correction) (output: RES1.DLNMselect.i.exp3.backward...)")


for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES1.DLNMselect.i.exp3[[i]])/5)
  RES1.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3[[i]]$expo_name), RES1.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES1.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES1.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp3[[i]]$var))),]
  RES1.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3[[i]]$CI.inf, RES1.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.bonferroni<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  #  RES1.DLNMselect.i.exp3[[i]]<-RES1.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.v2_backward<-RES1.DLNMselect.i.exp3

#save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

#*#
#load(paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES1.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp3.v2_backward[[i]]$p.value<0.05,1,0)
 
}

save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),"_FINAL.RData"))
#*#

print("Step 14 - Completed model penalized DLNM with variable selection Backward prunning with additional columns  (BH correction)  (output: RES1.DLNMselect.i.exp3.backward.v2...)")




###########################################################
###### MODEL 7: BACKWARD BY CORRECTION  	       ######
###########################################################


####  Penalized DLNM with variable selection and backward prunning (Xavier's proposal)

padj <- "by" #This argument allows to make more difficult a predictor to stay in the model during the backward prunning (correcting pvalues by the number of predictors in the model)

RES1.DLNMselect.i.exp3 <- vector("list", nsim)

RES1.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
  "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
  
  data.frame(applyDLNMbackwardsel(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,padj="by"),i)[,-10]
  
}

#save(RES1.DLNMselect.i.exp3,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.",as.character(padj),"_FINAL.RData"))

print("Step 15 - Completed model penalized DLNM with variable selection Backward prunning (BY correction) (output: RES1.DLNMselect.i.exp3.backward...)")


for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES1.DLNMselect.i.exp3[[i]])/5)
  RES1.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3[[i]]$expo_name), RES1.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES1.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES1.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES1.DLNMselect.i.exp3[[i]] <- RES1.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES1.DLNMselect.i.exp3[[i]]$var))),]
  RES1.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES1.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3[[i]]$CI.inf, RES1.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(rep(list.pval,each=5)<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.bonferroni<0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BH <0.05,1,0))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,ifelse(list.pval.BY <0.05,1,0))
  
  #  RES1.DLNMselect.i.exp3[[i]]<-RES1.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.v2_backward <- RES1.DLNMselect.i.exp3

#save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

#*#
#load(paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES1.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic <- ifelse(RES1.DLNMselect.i.exp3.v2_backward[[i]]$p.value<0.05,1,0)
 
}

save(RES1.DLNMselect.i.exp3.v2_backward,file=paste0("dataY1andX/exp3/RES1.DLNMselect.i.exp3.backward.v2.",as.character(padj),"_FINAL.RData"))
#*#

print("Step 16 - Completed model penalized DLNM with variable selection Backward prunning with additional columns (BY correction) (output: RES1.DLNMselect.i.exp3.backward.v2...)")



