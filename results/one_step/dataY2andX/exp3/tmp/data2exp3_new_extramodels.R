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
here::i_am("dataY2andX/exp3/data2exp3_new.R") # Declare the location of the current script and set wd / The project root is initialized (one step up), consistent with the location of the current script or report
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
load(file="dataY2andX/exp3/resu.sim.data2.exp3.RData")

# Load ExWAS selected variables for initialization
load(file="AVG_ExWAS_results/RES2av.Ewas.i.exp3.100.RData")

print("Step 2 - Loaded simulated predictor and outcome datasets")


################################# RUN MODELS

#### Setting up the parallelization process

parallel::detectCores()
n.cores <- 24#parallel::detectCores() - 1
my.cluster <- parallel::makeCluster( # create the cluster
  n.cores, 
  type = "PSOCK"
)
print(my.cluster) # check cluster definition (optional)
doParallel::registerDoParallel(cl = my.cluster) # register it to be used by %dopar%
foreach::getDoParRegistered() # check if it is registered (optional)
foreach::getDoParWorkers() # how many workers are available? (optional)

print("Step 3 - Set up for parallelization completed")


######## MODEL STARTED WITH VARIABLES SELECTED IN ExWAS Averaged APPROACH:
####  Penalized DLNM with variable selection and backward prunning (Xavier's proposal)

# Filtering input variables to those selected by AVG ExWAS

Selected_Xs <- vector("list", nsim)
for (i in 1:nsim) {
  
init_varlist <- gsub("_mean","_",RES2av.Ewas.i.exp3.100[[i]]$var[which(RES2av.Ewas.i.exp3.100[[i]]$EWAS.TP.none==1)])
init_varlist 
 ### eliminar todo tras .

Selected_Xs[[i]] <- resu.sim.dataX.i[[i]]$X[,unlist(lapply(as.character(init_varlist),function(x){grep(x, gsub("\\..*","_",colnames(resu.sim.dataX.i[[i]]$X)) )}))]

}

padj <- "none" #This argument allows to make more difficult a predictor to stay in the model during the backward prunning (correcting pvalues by the number of predictors in the model)

RES2.DLNMselect.i.exp3 <- vector("list", nsim)

RES2.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
                                                         "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
                                                           
                                                           data.frame(applyDLNMbackwardselAVG(data.Y=resu.sim.data2.exp3[[i]],data.X=Selected_Xs[[i]],padj="none"),i)[,-10]
                                                           
                                                         }

save(RES2.DLNMselect.i.exp3,file=paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.backward.AVGEXWAS.",as.character(padj),".RData"))

print("Step 4 - Completed model penalized DLNM with variable selection Backward prunning (output: RES2.DLNMselect.i.exp3.backward.AVGEXWAS...)")

for(i in 1:nsim) {
  
  RES2.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES2.DLNMselect.i.exp3[[i]])/5)
  RES2.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES2.DLNMselect.i.exp3[[i]]$expo_name), RES2.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES2.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES2.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES2.DLNMselect.i.exp3[[i]] <- RES2.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES2.DLNMselect.i.exp3[[i]] <- RES2.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES2.DLNMselect.i.exp3[[i]]$var))),]
  RES2.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES2.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data2.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES2.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES2.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES2.DLNMselect.i.exp3[[i]]$CI.inf, RES2.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                      ifelse(rep(list.pval,each=5)<0.05,1,0))
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                      ifelse(list.pval.bonferroni<0.05,1,0))
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                    ifelse(list.pval.BH <0.05,1,0))
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                    ifelse(list.pval.BY <0.05,1,0))
  #  RES2.DLNMselect.i.exp3[[i]]<-RES2.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES2.DLNMselect.i.exp3.v2_backward<-RES2.DLNMselect.i.exp3

save(RES2.DLNMselect.i.exp3.v2_backward,file=paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.backward.AVGEXWAS.v2.",as.character(padj),".RData"))

#*#
load(paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.backward.AVGEXWAS.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES2.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic <- ifelse(RES2.DLNMselect.i.exp3.v2_backward[[i]]$p.value<0.05,1,0)
  #RES2.DLNMselect.i.exp3.v2_backward[[i]]$DLNMpen.TP_dic[is.na(RES2.DLNMselect.i.exp3.v2_backward[[i]]$p.value)] <- 0
}

save(RES2.DLNMselect.i.exp3.v2_backward,file=paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.backward.AVGEXWAS.v2.",as.character(padj),".RData"))
#*#

print("Step 5 - Completed model penalized DLNM with variable selection Backward prunning with additional columns (output: RES2.DLNMselect.i.exp3.backward.AVGEXWAS.v2...)")








######## Penalized DLNM with Forward Stepwise Variable Selection FORWARD

padj <- "none" #This argument allows to make more difficult a predictor to stay in the model during the Forward prunning (correcting pvalues by the number of predictors in the model)

RES2.DLNMselect.i.exp3 <- vector("list", nsim)

RES2.DLNMselect.i.exp3 <- foreach(i = 1:nsim,.packages=c("foreach","Matrix","MASS","parallel","glmnet","spls","MXM","dlnm","splines","mgcv","doParallel",
                                                         "ranger","palmerpenguins","tidyverse","kableExtra")) %dopar% {
                                                           
                                                           data.frame(applyDLNMforwardsel(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,padj="none",i),i)[,-10]
                                                           
                                                         }

save(RES2.DLNMselect.i.exp3,file=paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.Forward.",as.character(padj),".RData"))

print("Step 6 - Completed model penalized DLNM with variable selection Forward prunning (output: RES2.DLNMselect.i.exp3.Forward...)")


for(i in 1:nsim) {
  
  RES2.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,nrow(RES2.DLNMselect.i.exp3[[i]])/5)
  RES2.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES2.DLNMselect.i.exp3[[i]]$expo_name), RES2.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  whole_dataset <- as.data.frame(names(resu.sim.dataX.i[[i]]$X))
  names(whole_dataset) <- "var"
  RES2.DLNMselect.i.exp3[[i]] <- merge.data.frame(whole_dataset,RES2.DLNMselect.i.exp3[[i]],by="var",all.x = TRUE)
  RES2.DLNMselect.i.exp3[[i]] <- RES2.DLNMselect.i.exp3[[i]][,c(2:9,1,10:11)]
  RES2.DLNMselect.i.exp3[[i]] <- RES2.DLNMselect.i.exp3[[i]][order(as.numeric(gsub("X","",RES2.DLNMselect.i.exp3[[i]]$var))),]
  RES2.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES2.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data2.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES2.DLNMselect.i.exp3[[i]]$p.value[grep(".1",RES2.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.bonferroni<-rep(p.adjust(list.pval,"bonferroni"),each=5)
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES2.DLNMselect.i.exp3[[i]]$CI.inf, RES2.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                      ifelse(rep(list.pval,each=5)<0.05,1,0))
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                      ifelse(list.pval.bonferroni<0.05,1,0))
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                    ifelse(list.pval.BH <0.05,1,0))
  RES2.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES2.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                    ifelse(list.pval.BY <0.05,1,0))
  #  RES2.DLNMselect.i.exp3[[i]]<-RES2.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES2.DLNMselect.i.exp3.v2_Forward<-RES2.DLNMselect.i.exp3

save(RES2.DLNMselect.i.exp3.v2_Forward,file=paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.Forward.v2.",as.character(padj),".RData"))

#*#
load(paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.Forward.v2.",as.character(padj),".RData"))

for(i in 1:nsim) {
  RES2.DLNMselect.i.exp3.v2_Forward[[i]]$DLNMpen.TP_dic <- ifelse(RES2.DLNMselect.i.exp3.v2_Forward[[i]]$p.value<0.05,1,0)
  #RES2.DLNMselect.i.exp3.v2_Forward[[i]]$DLNMpen.TP_dic[is.na(RES2.DLNMselect.i.exp3.v2_Forward[[i]]$p.value)] <- 0
}

save(RES2.DLNMselect.i.exp3.v2_Forward,file=paste0("dataY2andX/exp3/RES2.DLNMselect.i.exp3.Forward.v2.",as.character(padj),".RData"))
#*#

print("Step 7 - Completed model penalized DLNM with variable selection Forward prunning with additional columns (output: RES2.DLNMselect.i.exp3.Forward.v2...)")




