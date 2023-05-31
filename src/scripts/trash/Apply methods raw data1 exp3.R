library(Matrix) 
library(MASS) 
library(parallel) 
library(glmnet)
library(spls) 
library(MXM)
library(DSA)
library(dlnm)
library(splines)
library(mgcv)

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp3")
source(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/baseline/Functions method simulation.R")

# define number of simulated dataset (X and Y) and number of true X predictor
nsim<-100
nexp<-3

# load X data
load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/baseline/resu.sim.dataX.i.RData")

# load Y data
load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp3/resu.sim.data1.exp3.RData")


################################# RUN MODELS


##### EWAS

RES1.Ewas.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.Ewas.i.exp3[[i]] = applyEwas(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES1.Ewas.i.exp3[[i]]$numsim<-i
}

save(RES1.Ewas.i.exp3,file="RES1.Ewas.i.exp3.RData")
load(file="RES1.Ewas.i.exp3.RData")

##### ENET

init<-RES1.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.Enet.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.Enet.i.exp3[[i]] = applyEnet(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES1.Enet.i.exp3[[i]]$numsim<-i
}

save(RES1.Enet.i.exp3,file="RES1.Enet.i.exp3.RData")


##### sPLS

init<-RES1.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.sPLS.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.sPLS.i.exp3[[i]] = applySPLS(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES1.sPLS.i.exp3[[i]]$numsim<-i
}

save(RES1.sPLS.i.exp3,file="RES1.sPLS.i.exp3.RData")

##### MMPC

init<-RES1.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.MMPC.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.MMPC.i.exp3[[i]] = applyMMPC(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X)
  RES1.MMPC.i.exp3[[i]]$numsim<-i
}

save(RES1.MMPC.i.exp3,file="RES1.MMPC.i.exp3.RData")

##### DSA

init<-RES1.Ewas.i.exp3[[1]][,c("var","numsim")]  # store in init exposure names and true estimate
RES1.DSA.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.DSA.i.exp3[[i]] = applyDSA(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]]$X,maxsize=25)
  RES1.DSA.i.exp3[[i]]$numsim<-i
}

save(RES1.DSA.i.exp3,file="RES1.DSA.i.exp3.RData")


#### Penalized DLM

# check the significance of the overall p-value (crossbasis) and the confidence interval of each lag : if overall p-value is significant 
# (with or without correction for multiple testing) and the CI doesn't include 0, then considered as "selected" otherwise not (denominator=500).

RES1.DLNMpen.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.DLNMpen.i.exp3[[i]] = applyDLNMpen(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]])
  RES1.DLNMpen.i.exp3[[i]]$numsim<-i
  RES1.DLNMpen.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMpen.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  # p-value of the crossbasis (overall p-value)
  RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP<-as.numeric(as.character(RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP))
  # confidence interval of the lag includes 0 or not
  RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP<-as.numeric(as.character(ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP=="TRUE",1,0)))
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP[grep(".1",RES1.DLNMpen.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  RES1.DLNMpen.i.exp3[[i]]$DLNMpen.none<-ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP %in% c(0),0,
                                                ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP<0.05,1,0))
  RES1.DLNMpen.i.exp3[[i]]$DLNMpen.bonf<-ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP %in% c(0),0,
                                                ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP<0.0005,1,0))
  RES1.DLNMpen.i.exp3[[i]]$DLNMpen.bh<-ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP %in% c(0),0,
                                              ifelse(list.pval.BH <0.05,1,0))
  RES1.DLNMpen.i.exp3[[i]]$DLNMpen.by<-ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP %in% c(0),0,
                                              ifelse(list.pval.BY <0.05,1,0))
  RES1.DLNMpen.i.exp3[[i]]<-RES1.DLNMpen.i.exp3[[i]][,c("var","numsim","true.pred","DLNMpen.TP","estimate","CI.inf","CI.sup","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by","DLNMpen_lag.TP")]
  
  # RES1.DLNMpen.i.exp3[[i]]$sign.lag<-ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP>0.05,0,
  #                                      ifelse(RES1.DLNMpen.i.exp3[[i]]$DLNMpen.TP<0.05 & RES1.DLNMpen.i.exp3[[i]]$DLNMpen_lag.TP==0,0,1))
}

save(RES1.DLNMpen.i.exp3,file="RES1.DLNMpen.i.exp3.RData")

 
#### Penalized DLM with variable selection

RES1.DLNMselect.i.exp3<- vector("list", nsim)
for(i in 1:nsim) {
  print(i)
  RES1.DLNMselect.i.exp3[[i]] = applyDLNMselect(data.Y=resu.sim.data1.exp3[[i]],data.X=resu.sim.dataX.i[[i]])
  RES1.DLNMselect.i.exp3[[i]]$numsim<-i
}

save(RES1.DLNMselect.i.exp3,file="RES1.DLNMselect.i.exp3.RData")


for(i in 1:nsim) {
  
  RES1.DLNMselect.i.exp3[[i]]$num_time<-rep(1:5,100)
  RES1.DLNMselect.i.exp3[[i]]$var<-paste(toupper(RES1.DLNMselect.i.exp3[[i]]$expo_name), RES1.DLNMselect.i.exp3[[i]]$num_time,sep = ".")
  
  RES1.DLNMselect.i.exp3[[i]]$true.pred<-ifelse(RES1.DLNMselect.i.exp3[[i]]$var %in% resu.sim.data1.exp3[[i]]$true.pred,1,0)
  
  # store p-values (by step of 5 = overall p-value of the crossbasis)
  list.pval<-RES1.DLNMselect.i.exp3[[i]]$`p-value`[grep(".1",RES1.DLNMselect.i.exp3[[i]]$var,fixed = TRUE)]
  list.pval.BH<-rep(p.adjust(list.pval,"BH"),each=5)
  list.pval.BY<-rep(p.adjust(list.pval,"BY"),each=5)
  
  # CI includes 0 or not
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-!(data.table::between(0, RES1.DLNMselect.i.exp3[[i]]$CI.inf, RES1.DLNMselect.i.exp3[[i]]$CI.sup))
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP==FALSE,0,1)
  
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.none<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                      ifelse(list.pval<0.05,1,0))
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bonf<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                      ifelse(list.pval<0.0005,1,0))
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.bh<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                    ifelse(list.pval.BH <0.05,1,0))
  RES1.DLNMselect.i.exp3[[i]]$DLNMselect.by<-ifelse(RES1.DLNMselect.i.exp3[[i]]$DLNMselect_lag.TP %in% c(0),0,
                                                    ifelse(list.pval.BY <0.05,1,0))
  #  RES1.DLNMselect.i.exp3[[i]]<-RES1.DLNMselect.i.exp3[[i]][,c("var","numsim","true.pred","DLNMselect.TP","estimate","CI.inf","CI.sup","DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by","DLNMselect_lag.TP")]
  
}

RES1.DLNMselect.i.exp3.v2<-RES1.DLNMselect.i.exp3

save(RES1.DLNMselect.i.exp3.v2,file="RES1.DLNMselect.i.exp3.v2.RData")



