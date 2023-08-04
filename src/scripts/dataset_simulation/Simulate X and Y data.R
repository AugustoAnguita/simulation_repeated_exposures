library(MASS)
library(Matrix)
library(dplyr)
library(matrixcalc)


###### Script to simulate the data under the different scenarios

# load the functions
source(file="./src/source_functions/functions.R")

# define number of simulated dataset (X and Y)
nsim<-100


###### Simulate X data

# load the correlation matrix to be used to simulate X
load(file="./data/Xs/Sigma.Rdata") # Real correlation matrix observed in Helix project using postnatal exposure data

# Simulate X data based on the correlation matrix - 100 times 
resu.sim.dataX.i <- vector("list", nsim)
for(i in 1:nsim) {
   resu.sim.dataX.i[[i]] = simX(Sigma=Sigma, N = 1200)
}

# reordering columns
col_order<-c("X1.1","X1.2","X1.3","X1.4","X1.5","X2.1","X2.2","X2.3","X2.4","X2.5","X3.1","X3.2","X3.3","X3.4","X3.5","X4.1","X4.2","X4.3","X4.4","X4.5","X5.1","X5.2","X5.3","X5.4","X5.5","X6.1","X6.2","X6.3","X6.4","X6.5","X7.1","X7.2","X7.3","X7.4","X7.5","X8.1","X8.2","X8.3","X8.4","X8.5","X9.1","X9.2","X9.3","X9.4","X9.5",
             "X10.1","X10.2","X10.3","X10.4","X10.5","X11.1","X11.2","X11.3","X11.4","X11.5","X12.1","X12.2","X12.3","X12.4","X12.5","X13.1","X13.2","X13.3","X13.4","X13.5","X14.1","X14.2","X14.3","X14.4","X14.5","X15.1","X15.2","X15.3","X15.4","X15.5","X16.1","X16.2","X16.3","X16.4","X16.5","X17.1","X17.2","X17.3","X17.4","X17.5","X18.1","X18.2","X18.3","X18.4","X18.5","X19.1","X19.2","X19.3","X19.4","X19.5",
             "X20.1","X20.2","X20.3","X20.4","X20.5","X21.1","X21.2","X21.3","X21.4","X21.5","X22.1","X22.2","X22.3","X22.4","X22.5","X23.1","X23.2","X23.3","X23.4","X23.5","X24.1","X24.2","X24.3","X24.4","X24.5","X25.1","X25.2","X25.3","X25.4","X25.5","X26.1","X26.2","X26.3","X26.4","X26.5","X27.1","X27.2","X27.3","X27.4","X27.5","X28.1","X28.2","X28.3","X28.4","X28.5","X29.1","X29.2","X29.3","X29.4","X29.5",
             "X30.1","X30.2","X30.3","X30.4","X30.5","X31.1","X31.2","X31.3","X31.4","X31.5","X32.1","X32.2","X32.3","X32.4","X32.5","X33.1","X33.2","X33.3","X33.4","X33.5","X34.1","X34.2","X34.3","X34.4","X34.5","X35.1","X35.2","X35.3","X35.4","X35.5","X36.1","X36.2","X36.3","X36.4","X36.5","X37.1","X37.2","X37.3","X37.4","X37.5","X38.1","X38.2","X38.3","X38.4","X38.5","X39.1","X39.2","X39.3","X39.4","X39.5",
             "X40.1","X40.2","X40.3","X40.4","X40.5","X41.1","X41.2","X41.3","X41.4","X41.5","X42.1","X42.2","X42.3","X42.4","X42.5","X43.1","X43.2","X43.3","X43.4","X43.5","X44.1","X44.2","X44.3","X44.4","X44.5","X45.1","X45.2","X45.3","X45.4","X45.5","X46.1","X46.2","X46.3","X46.4","X46.5","X47.1","X47.2","X47.3","X47.4","X47.5","X48.1","X48.2","X48.3","X48.4","X48.5","X49.1","X49.2","X49.3","X49.4","X49.5",
             "X50.1","X50.2","X50.3","X50.4","X50.5","X51.1","X51.2","X51.3","X51.4","X51.5","X52.1","X52.2","X52.3","X52.4","X52.5","X53.1","X53.2","X53.3","X53.4","X53.5","X54.1","X54.2","X54.3","X54.4","X54.5","X55.1","X55.2","X55.3","X55.4","X55.5","X56.1","X56.2","X56.3","X56.4","X56.5","X57.1","X57.2","X57.3","X57.4","X57.5","X58.1","X58.2","X58.3","X58.4","X58.5","X59.1","X59.2","X59.3","X59.4","X59.5",
             "X60.1","X60.2","X60.3","X60.4","X60.5","X61.1","X61.2","X61.3","X61.4","X61.5","X62.1","X62.2","X62.3","X62.4","X62.5","X63.1","X63.2","X63.3","X63.4","X63.5","X64.1","X64.2","X64.3","X64.4","X64.5","X65.1","X65.2","X65.3","X65.4","X65.5","X66.1","X66.2","X66.3","X66.4","X66.5","X67.1","X67.2","X67.3","X67.4","X67.5","X68.1","X68.2","X68.3","X68.4","X68.5","X69.1","X69.2","X69.3","X69.4","X69.5",
             "X70.1","X70.2","X70.3","X70.4","X70.5","X71.1","X71.2","X71.3","X71.4","X71.5","X72.1","X72.2","X72.3","X72.4","X72.5","X73.1","X73.2","X73.3","X73.4","X73.5","X74.1","X74.2","X74.3","X74.4","X74.5","X75.1","X75.2","X75.3","X75.4","X75.5","X76.1","X76.2","X76.3","X76.4","X76.5","X77.1","X77.2","X77.3","X77.4","X77.5","X78.1","X78.2","X78.3","X78.4","X78.5","X79.1","X79.2","X79.3","X79.4","X79.5",
             "X80.1","X80.2","X80.3","X80.4","X80.5","X81.1","X81.2","X81.3","X81.4","X81.5","X82.1","X82.2","X82.3","X82.4","X82.5","X83.1","X83.2","X83.3","X83.4","X83.5","X84.1","X84.2","X84.3","X84.4","X84.5","X85.1","X85.2","X85.3","X85.4","X85.5","X86.1","X86.2","X86.3","X86.4","X86.5","X87.1","X87.2","X87.3","X87.4","X87.5","X88.1","X88.2","X88.3","X88.4","X88.5","X89.1","X89.2","X89.3","X89.4","X89.5",
             "X90.1","X90.2","X90.3","X90.4","X90.5","X91.1","X91.2","X91.3","X91.4","X91.5","X92.1","X92.2","X92.3","X92.4","X92.5","X93.1","X93.2","X93.3","X93.4","X93.5","X94.1","X94.2","X94.3","X94.4","X94.5","X95.1","X95.2","X95.3","X95.4","X95.5","X96.1","X96.2","X96.3","X96.4","X96.5","X97.1","X97.2","X97.3","X97.4","X97.5","X98.1","X98.2","X98.3","X98.4","X98.5","X99.1","X99.2","X99.3","X99.4","X99.5",
             "X100.1","X100.2","X100.3","X100.4","X100.5")

for(i in 1:nsim) {
  resu.sim.dataX.i[[i]]$X = resu.sim.dataX.i[[i]]$X[,col_order]
}

# Save X data
save(resu.sim.dataX.i,file="./data/Xs/resu.sim.dataX.i.RData")



###### Simulate Y data under scenario 1

# Simulate Y data (simtestY1) : all X time points are causaly associated with Y
#     X = Exposure data
#     nsim = number of simulations = 100
#     n.true.exposures = k = number of exposures associated with Y (multiplied by number of time points) = 3*5
#     nvars = Number of exposure for each time point = 100
#     N = number of subjects = 1200


### Simulate Y data, under scenario 1, when k=3 
nexp<-3

resu.sim.data1.exp3 <- vector("list", nsim)
for(i in 1:nsim) {
  data.X<-as.data.frame(resu.sim.dataX.i[[i]]$X)
  resu.sim.data1.exp3[[i]] = simtestY1(X=data.X, n.true.exposures=nexp, nvars=100, N = 1200)
}

# Save Y data, under scenario 1, when k=3
save(resu.sim.data1.exp3,file="./data/Ys/data1exp3/resu.sim.data1.exp3.RData")


### Simulate Y data, under scenario 1, when k=5 
nexp<-5

resu.sim.data1.exp5 <- vector("list", nsim)
for(i in 1:nsim) {
  data.X<-as.data.frame(resu.sim.dataX.i[[i]]$X)
  resu.sim.data1.exp5[[i]] = simtestY1(X=data.X,n.true.exposures=nexp,nvars=100,N = 1200)
}

# Save Y data, under scenario 1, when k=5
save(resu.sim.data1.exp5,file="./data/Ys/data1exp5/resu.sim.data1.exp5.RData")


### Simulate Y data, under scenario 1, when k=10
nexp<-10

resu.sim.data1.exp10 <- vector("list", nsim)
for(i in 1:nsim) {
  data.X<-as.data.frame(resu.sim.dataX.i[[i]]$X)
  resu.sim.data1.exp10[[i]] = simtestY1(X=data.X,n.true.exposures=nexp,nvars=100,N = 1200)
}

# Save Y data, under scenario 1, when k=10
save(resu.sim.data1.exp10,file="./data/Ys/data1exp10/resu.sim.data1.exp10.RData")


###### Simulate Y data under scenario 2

# Simulate Y data (simtestY2) : X measured at all time point associated with Y
#     X = Exposure data
#     nsim = number of simulation = 100
#     n.true.exposures = k = number of exposures associated with Y (multiplied by number of time points = 1) = 3*1
#     nvars = Number of exposure for each time point = 100
#     N = number of subjects = 1200



### Simulate Y data, under scenario 2, when k=3 
nexp<-3

resu.sim.data2.exp3 <- vector("list", nsim)
for(i in 1:nsim) {
  data.X<-as.data.frame(resu.sim.dataX.i[[i]]$X)
  resu.sim.data2.exp3[[i]] = simtestY2(X=data.X,n.true.exposures=nexp,nvars=100,N = 1200)
}

# Save Y data, under scenario 2, when k=3
save(resu.sim.data2.exp3,file="./data/Ys/data2exp3/resu.sim.data2.exp3.RData")



### Simulate Y data, under scenario 2, when k=5 
nexp<-5

resu.sim.data2.exp5 <- vector("list", nsim)
for(i in 1:nsim) {
  data.X<-as.data.frame(resu.sim.dataX.i[[i]]$X)
  resu.sim.data2.exp5[[i]] = simtestY2(X=data.X,n.true.exposures=nexp,nvars=100,N = 1200)
}

# Save Y data, under scenario 2, when k=5
save(resu.sim.data2.exp5,file="./data/Ys/data2exp3/resu.sim.data2.exp5.RData")


### Simulate Y data, under scenario 2, when k=10
nexp<-10

resu.sim.data2.exp10 <- vector("list", nsim)
for(i in 1:nsim) {
  data.X<-as.data.frame(resu.sim.dataX.i[[i]]$X)
  resu.sim.data2.exp10[[i]] = simtestY2(X=data.X,n.true.exposures=nexp,nvars=100,N = 1200)
}

# Save Y data, under scenario 2, when k=10
save(resu.sim.data2.exp10,file="./data/Ys/data2exp3/resu.sim.data2.exp10.RData")






################## Reorganizing X data for the 2-step approaches



##### Reshape dataX from wide to long format

reshap.long<-function(X,nvars,r){
  lnames = list()
  for (j in 1:nvars) {
    lnames[[j]] = matrix(names(X),ncol=r,byrow=T)[j,] 
  }
  
  Xlong = reshape(X,direction="long",varying=lnames)
  Xlong = Xlong[,c(dim(Xlong)[2],1:(dim(Xlong)[2]-1))]
  names(Xlong) = c("id","time",paste0("X",1:nvars)) 
  # Xlong$id<-as.character(Xlong$id)
  Xlong = Xlong[order(Xlong$id,Xlong$time),]
  return(Xlong)
}

Xlong.i <- vector("list", nsim)
for(i in 1:nsim) {
  Xlong.i[[i]]<-reshap.long(X=resu.sim.dataX.i[[i]]$X,nvars=100,r=5)
}

save(Xlong.i,file="./data/Xs/Xlong.i.RData")


##### Calculate the average exposure level across the time points

average.long<-function(data.X){
  data.X.av<- data.X %>%
    group_by((id)) %>%
    summarise_at(.vars = c(3:102),
                 .funs = c(mean="mean"))
  data.X.av<-data.X.av[-1]
  return(data.X.av)
}

Xave.i <- vector("list", nsim)
for(i in 1:nsim) {
  Xave.i[[i]]<-average.long(data.X=Xlong.i[[i]])
}

save(Xave.i,file="./data/Xs/Xave.i.RData")


##### Renaming the true predictors of Y stored in resu.sim.dataX.expX$true.pred (from .1, .2, .3, .4, .5 to _mean)

### Under scenario 1, when k=3
nexp<-3

# load Y data (simulated from raw X data, not X average)
load(file="./data/Ys/data1exp3/resu.sim.data1.exp3.RData")

# create a new object with Y data modifying the true.pred list (eg, from X15.1 X15.2 X15.3 X15.4 X15.5 to X15_mean)
resu.sim.data1av.exp3.100<-resu.sim.data1.exp3

for (i in 1:100){
  resu.sim.data1av.exp3.100[[i]]$true.pred<-paste0(ifelse(nchar(as.character(resu.sim.data1av.exp3.100[[i]]$true.pred))==4,substring(resu.sim.data1av.exp3.100[[i]]$true.pred,first=1,last=2),
                                                      ifelse(nchar(as.character(resu.sim.data1av.exp3.100[[i]]$true.pred))==5,substring(resu.sim.data1av.exp3.100[[i]]$true.pred,first=1,last=3),
                                                             substring(resu.sim.data1av.exp3.100[[i]]$true.pred,first=1,last=4))),"_mean")
  resu.sim.data1av.exp3.100[[i]]$true.pred<-resu.sim.data1av.exp3.100[[i]]$true.pred[seq(1,5*nexp,by=5)]
  resu.sim.data1av.exp3.100[[i]]$beta<-resu.sim.data1av.exp3.100[[i]]$beta[seq(1,5*nexp,by=5)]
  
}

save(resu.sim.data1av.exp3.100,file="./data/Ys/data1exp3/twostep/resu.sim.data1av.exp3.100.RData")


### Under scenario 1, when k=5
nexp<-5

# load Y data (simulated from raw X data, not X average)
load(file="./data/Ys/data1exp5/resu.sim.data1.exp5.RData")

# create a new object with Y data modifying the true.pred list (eg, from X15.1 X15.2 X15.3 X15.4 X15.5 to X15_mean)
resu.sim.data1av.exp5.100<-resu.sim.data1.exp5

for (i in 1:100){
  resu.sim.data1av.exp5.100[[i]]$true.pred<-paste0(ifelse(nchar(as.character(resu.sim.data1av.exp5.100[[i]]$true.pred))==4,substring(resu.sim.data1av.exp5.100[[i]]$true.pred,first=1,last=2),
                                                          ifelse(nchar(as.character(resu.sim.data1av.exp5.100[[i]]$true.pred))==5,substring(resu.sim.data1av.exp5.100[[i]]$true.pred,first=1,last=3),
                                                                 substring(resu.sim.data1av.exp5.100[[i]]$true.pred,first=1,last=4))),"_mean")
  resu.sim.data1av.exp5.100[[i]]$true.pred<-resu.sim.data1av.exp5.100[[i]]$true.pred[seq(1,5*nexp,by=5)]
  resu.sim.data1av.exp5.100[[i]]$beta<-resu.sim.data1av.exp5.100[[i]]$beta[seq(1,5*nexp,by=5)]
  
}

save(resu.sim.data1av.exp5.100,file="./data/Ys/data1exp5/twostep/resu.sim.data1av.exp5.100.RData")


### Under scenario 1, when k=10
nexp<-10

# load Y data (simulated from raw X data, not X average)
load(file="./data/Ys/data1exp10/resu.sim.data1.exp10.RData")

# create a new object with Y data modifying the true.pred list (eg, from X15.1 X15.2 X15.3 X15.4 X15.5 to X15_mean)
resu.sim.data1av.exp10.100<-resu.sim.data1.exp10

for (i in 1:100){
  resu.sim.data1av.exp10.100[[i]]$true.pred<-paste0(ifelse(nchar(as.character(resu.sim.data1av.exp10.100[[i]]$true.pred))==4,substring(resu.sim.data1av.exp10.100[[i]]$true.pred,first=1,last=2),
                                                          ifelse(nchar(as.character(resu.sim.data1av.exp10.100[[i]]$true.pred))==5,substring(resu.sim.data1av.exp10.100[[i]]$true.pred,first=1,last=3),
                                                                 substring(resu.sim.data1av.exp10.100[[i]]$true.pred,first=1,last=4))),"_mean")
  resu.sim.data1av.exp10.100[[i]]$true.pred<-resu.sim.data1av.exp10.100[[i]]$true.pred[seq(1,5*nexp,by=5)]
  resu.sim.data1av.exp10.100[[i]]$beta<-resu.sim.data1av.exp10.100[[i]]$beta[seq(1,5*nexp,by=5)]
  
}

save(resu.sim.data1av.exp10.100,file="./data/Ys/data1exp10/twostep/resu.sim.data1av.exp10.100.RData")


### Under scenario 2, when k=3
nexp<-3

# load Y data (simulated from raw X data, not X average)
load(file="./data/Ys/data2exp3/resu.sim.data2.exp3.RData")

# create a new object with Y data modifying the true.pred list (eg, from X15.1 X15.2 X15.3 X15.4 X15.5 to X15_mean)
resu.sim.data2av.exp3.100<-resu.sim.data2.exp3

for (i in 1:100){
  resu.sim.data2av.exp3.100[[i]]$true.pred<-paste0(ifelse(nchar(as.character(resu.sim.data2av.exp3.100[[i]]$true.pred))==4,substring(resu.sim.data2av.exp3.100[[i]]$true.pred,first=1,last=2),
                                                          ifelse(nchar(as.character(resu.sim.data2av.exp3.100[[i]]$true.pred))==5,substring(resu.sim.data2av.exp3.100[[i]]$true.pred,first=1,last=3),
                                                                 substring(resu.sim.data2av.exp3.100[[i]]$true.pred,first=1,last=4))),"_mean")
  resu.sim.data2av.exp3.100[[i]]$true.pred<-resu.sim.data2av.exp3.100[[i]]$true.pred[seq(1,5*nexp,by=5)]
  resu.sim.data2av.exp3.100[[i]]$beta<-resu.sim.data2av.exp3.100[[i]]$beta[seq(1,5*nexp,by=5)]
  
}

save(resu.sim.data2av.exp3.100,file="./data/Ys/data2exp3/twostep/resu.sim.data2av.exp3.100.RData")


### Under scenario 2, when k=5
nexp<-5

# load Y data (simulated from raw X data, not X average)
load(file="./data/Ys/data2exp5/resu.sim.data2.exp5.RData")

# create a new object with Y data modifying the true.pred list (eg, from X15.1 X15.2 X15.3 X15.4 X15.5 to X15_mean)
resu.sim.data2av.exp5.100<-resu.sim.data2.exp5

for (i in 1:100){
  resu.sim.data2av.exp5.100[[i]]$true.pred<-paste0(ifelse(nchar(as.character(resu.sim.data2av.exp5.100[[i]]$true.pred))==4,substring(resu.sim.data2av.exp5.100[[i]]$true.pred,first=1,last=2),
                                                          ifelse(nchar(as.character(resu.sim.data2av.exp5.100[[i]]$true.pred))==5,substring(resu.sim.data2av.exp5.100[[i]]$true.pred,first=1,last=3),
                                                                 substring(resu.sim.data2av.exp5.100[[i]]$true.pred,first=1,last=4))),"_mean")
  resu.sim.data2av.exp5.100[[i]]$true.pred<-resu.sim.data2av.exp5.100[[i]]$true.pred[seq(1,5*nexp,by=5)]
  resu.sim.data2av.exp5.100[[i]]$beta<-resu.sim.data2av.exp5.100[[i]]$beta[seq(1,5*nexp,by=5)]
  
}

save(resu.sim.data2av.exp5.100,file="./data/Ys/data2exp5/twostep/resu.sim.data2av.exp5.100.RData")


### Under scenario 2, when k=10
nexp<-10

# load Y data (simulated from raw X data, not X average)
load(file="./data/Ys/data2exp10/resu.sim.data2.exp10.RData")

# create a new object with Y data modifying the true.pred list (eg, from X15.1 X15.2 X15.3 X15.4 X15.5 to X15_mean)
resu.sim.data2av.exp10.100<-resu.sim.data2.exp10

for (i in 1:100){
  resu.sim.data2av.exp10.100[[i]]$true.pred<-paste0(ifelse(nchar(as.character(resu.sim.data2av.exp10.100[[i]]$true.pred))==4,substring(resu.sim.data2av.exp10.100[[i]]$true.pred,first=1,last=2),
                                                          ifelse(nchar(as.character(resu.sim.data2av.exp10.100[[i]]$true.pred))==5,substring(resu.sim.data2av.exp10.100[[i]]$true.pred,first=1,last=3),
                                                                 substring(resu.sim.data2av.exp10.100[[i]]$true.pred,first=1,last=4))),"_mean")
  resu.sim.data2av.exp10.100[[i]]$true.pred<-resu.sim.data2av.exp10.100[[i]]$true.pred[seq(1,5*nexp,by=5)]
  resu.sim.data2av.exp10.100[[i]]$beta<-resu.sim.data2av.exp10.100[[i]]$beta[seq(1,5*nexp,by=5)]
  
}

save(resu.sim.data2av.exp10.100,file="./data/Ys/data2exp10/twostep/resu.sim.data2av.exp10.100.RData")

