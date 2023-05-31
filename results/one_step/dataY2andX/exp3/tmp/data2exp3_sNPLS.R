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
here::i_am("dataY2andX/exp3/data2exp3_sNPLS.R") # Declare the location of the current script and set wd / The project root is initialized (one step up), consistent with the location of the current script or report
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



##### Apply sNPLS

library(future)
future::plan(cluster, workers=my.cluster)
RES2.sNPLS.i.exp3 <- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES2.sNPLS.i.exp3[[i]] = applysNPLS(data.Y=resu.sim.data2.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES2.sNPLS.i.exp3[[i]]$numsim <- i
}
save(RES2.sNPLS.i.exp3,file="RES2.sNPLS.i.exp3.RData")










################# MANUAL
#### sNPLS
# 
# colnames(resu.sim.dataX.i[[1]]$X)
# 
# Y <- resu.sim.data2.exp3[[1]]$resp
# X <- resu.sim.dataX.i[[1]]$X
# 
# library(stringr)
# names <- str_replace(colnames(X), '.+.(.+)', '\\1')
# 
# Ts <- c("1","2","3","4","5")
# 
# arr = array(NA, dim=c(nrow(X),ncol(X)/5,5))
# for (i in Ts) {
#   arr[,,as.numeric(i)] <- as.matrix(X[,grep(i,names)])
# }
# dim(arr)
# 
# library(sNPLS)
# library(future)
# future::plan(cluster, workers=my.cluster)
# cv_result <- repeat_cv(arr,Y,ncomp = 1:3, keepJ = 1:30, keepK = 1:5, parallel = TRUE, times = 5, nfold = 5)
# library(plyr)
# select_best_hp <- ddply(cv_result,.(ncomp,keepJ,keepK),nrow)
# if ( length(which(select_best_hp$V1 %in% max(select_best_hp$V1))) > 1 ) {
#   selected_hp <- select_best_hp[which(select_best_hp$V1 %in% max(select_best_hp$V1)),1:3]
#   selected_hp <- subset(selected_hp, ncomp > 1)
# }else{
# selected_hp <- select_best_hp[which(select_best_hp$V1 %in% max(select_best_hp$V1)),1:3]  
# }
# 
# fit <- sNPLS(arr,Y, ncomp = selected_hp[1,1], keepJ = rep(selected_hp[1,2],selected_hp[1,1]), keepK = rep(selected_hp[1,3],selected_hp[1,1]), silent = FALSE)
# 
# # Fig. 1 shows how the 1000 samples spread over the two ﬁrst compo-
# #   nents. Fig. 2 is similar, but related to the y scores.
# plot(fit, type="T", cex.axis=1.2, cex.lab=1.2, las=1, bty="L")
# plot(fit, type="U", cex.axis=1.2, cex.lab=1.2, las=1, bty="L")
# 
# # On the other hand, Fig. 3 Wj presents the attributes and Fig. 4 Wk the times
# # that help in predicting the y variable, for each of the two components.
# plot(fit, type="Wj", cex.axis=1.2, cex.lab=1.2, las=1, bty="L")
# plot(fit, type="Wk", cex.axis=1.2, cex.lab=1.2, las=1, bty="L")
# 
# # It can be seen how, for the ﬁrst component, there is approximately the
# # same effect of all times with respect variable 52 (related to
# # the ﬁrst component in the attributes mode) when trying to predict the
# # Y. In this case, since the times' weights show positive values,
# # the higher the value of attribute 52, the higher the Y. Among all times, it seems 
# # that time points 2 and 5 has a higher influence.
# # For the second component, we can observe the same, but
# # for attribute 32.
# plot(fit, type="variables", cex.axis=1.2, cex.lab=1.2, las=1, bty="L", lwd=2)
# plot(fit, type="time", cex.axis=1.2, cex.lab=1.2, las=1, bty="L", lwd=2, xlab="Time")
# 
# output_coefs <- as.data.frame(summary(fit))
# rownames(output_coefs) <- gsub("X.","X",rownames(output_coefs))
# names(output_coefs) <- c(1:5)
# output_coefs$var <- rownames(output_coefs)
# output_coefs <- output_coefs %>%
# pivot_longer(!var, names_to = "time", values_to = "coef")
# output_coefs$var <- paste(output_coefs$var,output_coefs$time, sep = ".")
# output_coefs <- output_coefs[,c(1,3)]
# output_coefs$sNPLS.TP <- as.factor(ifelse(is.na(output_coefs$coef),"0",ifelse(output_coefs$coef%in%"0","0","1")))
# output_coefs$var[which(output_coefs$sNPLS.TP==1)] == resu.sim.data2.exp3[[1]]$true.pred
# 
# #Check results:
# sort(as.character(output_coefs$var[which(output_coefs$sNPLS.TP==1)])) == sort(as.character(resu.sim.data2.exp3[[1]]$true.pred))




