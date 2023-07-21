

nsim<-100

####### recup ICC dataX 

load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/baseline/resu.sim.dataX.i.RData")

ICC.i<-vector("list",nsim)
for (i in 1:nsim){
  a<-NULL
  for (j in 1:100){
    a<-c(a,rep(resu.sim.dataX.i[[i]]$rho.vec[j],5))
  }
  var<-colnames(resu.sim.dataX.i[[i]]$X)
  numsim<-i
  ICC_var<-cbind(var,numsim,a)
  ICC.i[[i]]<-as.data.frame(rbind(ICC.i[[i]],ICC_var))
  #  ICC.i[[i]]<-ICC.i[[i]][-1,]
}

########################## Performance to identify the true exposure at the true time point ####################################


############################################# DATA1

############### EXP3

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp3")
nexp<-3


load(file="RES1.all.exp3.2023.RData")
colnames(RES1.all.exp3)
colnames(RES1.all.exp3)<-c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES1.all.exp3<-RES1.all.exp3[,c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                                "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]

load(file="twostep/RES1avRed.all.exp3.0621.RData")
colnames(RES1avRed.all.exp3.0621)
colnames(RES1avRed.all.exp3.0621)<-c("var","numsim","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES1.DLNM.all.exp3.RData")
colnames(RES1.DLNM.all.exp3)


load(file="dlnm/RES1.DLNM.AVG.all.exp3.RData")
colnames(RES1.DLNM.AVG.all.exp3)
colnames(RES1.DLNM.AVG.all.exp3)<-c("var","numsim","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                    "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES1.total.exp3<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                        list(RES1.all.exp3,RES1avRed.all.exp3.0621,RES1.DLNM.all.exp3,RES1.DLNM.AVG.all.exp3))
  
  
### Merge avec ICC

ICC.data1.exp3<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data1.exp3)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data1.exp3<-rbind(ICC.data1.exp3,ICC.i[[i]])
}
ICC.data1.exp3<-ICC.data1.exp3[-1,]
ICC.data1.exp3$numsim<-as.integer(ICC.data1.exp3$numsim)
colnames(ICC.data1.exp3)<-c("var","numsim","ICC")

RES1.total.exp3.ICC<-merge(RES1.total.exp3,ICC.data1.exp3,by=c("var","numsim"))

table(RES1.total.exp3.ICC$true.pred,RES1.total.exp3.ICC$ICC)



library(dplyr)

squareTable <- function(x,y) {
  x <- factor(x,levels = c("0","1"))
  y <- factor(y,levels = c("0","1"))
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}




RES1.total.exp3.ICC.low<-RES1.total.exp3.ICC[RES1.total.exp3.ICC$ICC==0.1,]
RES1.total.exp3.ICC.med<-RES1.total.exp3.ICC[RES1.total.exp3.ICC$ICC==0.5,]
RES1.total.exp3.ICC.high<-RES1.total.exp3.ICC[RES1.total.exp3.ICC$ICC==0.9,]


### All

RES1.i.total.exp3.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp3.ICC[[i]]<-RES1.total.exp3.ICC[RES1.total.exp3.ICC$numsim == i,]
}

result.data1.i.exp3<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp3.ICC[[i]])[2]-1)){
    RES1.i.total.exp3.ICC[[i]][,j]<-as.factor(as.character(RES1.i.total.exp3.ICC[[i]][,j]))  
    RES1.i.total.exp3.ICC[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp3.ICC[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp3.ICC[[i]][,j],RES1.i.total.exp3.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp3.ICC[[i]][,j],RES1.i.total.exp3.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp3.ICC[[i]][,j],RES1.i.total.exp3.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp3.ICC[[i]][,j],RES1.i.total.exp3.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp3.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp3[[i]]<-rbind(result.data1.i.exp3[[i]],tot)
   }
}


result.data1.exp3<-matrix(ncol=dim(result.data1.i.exp3[[1]])[2])
for (i in 1:nsim){
  result.data1.exp3<-rbind(result.data1.exp3,result.data1.i.exp3[[i]])
}

result.data1.exp3<-as.data.frame(result.data1.exp3)
colnames(result.data1.exp3)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3$nselected<-as.numeric(as.character(result.data1.exp3$nselected))
result.data1.exp3$tn<-as.numeric(as.character(result.data1.exp3$tn))
result.data1.exp3$fn<-as.numeric(as.character(result.data1.exp3$fn))
result.data1.exp3$fp<-as.numeric(as.character(result.data1.exp3$fp))
result.data1.exp3$tp<-as.numeric(as.character(result.data1.exp3$tp))

result.data1.exp3$sensitivity<-round(result.data1.exp3$tp/(result.data1.exp3$tp+result.data1.exp3$fn)*100,1)
result.data1.exp3$specificity<-round(result.data1.exp3$tn/(result.data1.exp3$tn+result.data1.exp3$fp)*100,1)
result.data1.exp3$false.pos.rate<-round(result.data1.exp3$fp/(result.data1.exp3$fp+result.data1.exp3$tn)*100,1)
result.data1.exp3$false.neg.rate<-round(result.data1.exp3$fn/(result.data1.exp3$fn+result.data1.exp3$tp)*100,1)
result.data1.exp3$false.disc.rate<-round(result.data1.exp3$fp/(result.data1.exp3$fp+result.data1.exp3$tp)*100,1)
result.data1.exp3$false.omit.rate<-round(result.data1.exp3$fn/(result.data1.exp3$fn+result.data1.exp3$tn)*100,1)

result.data1.exp3<-result.data1.exp3[-1,]

result.data1.exp3$expo<-3
result.data1.exp3$ICC<-1


# Low ICC

RES1.i.total.exp3.ICC.low<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp3.ICC.low[[i]]<-RES1.total.exp3.ICC.low[RES1.total.exp3.ICC.low$numsim == i,]
}

result.data1.i.exp3.ICC.low<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp3.ICC.low[[i]])[2]-1)){
    RES1.i.total.exp3.ICC.low[[i]][,j]<-as.factor(as.character(RES1.i.total.exp3.ICC.low[[i]][,j]))  
    RES1.i.total.exp3.ICC.low[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp3.ICC.low[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp3.ICC.low[[i]][,j],RES1.i.total.exp3.ICC.low[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp3.ICC.low[[i]][,j],RES1.i.total.exp3.ICC.low[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp3.ICC.low[[i]][,j],RES1.i.total.exp3.ICC.low[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp3.ICC.low[[i]][,j],RES1.i.total.exp3.ICC.low[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp3.ICC.low[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp3.ICC.low[[i]]<-rbind(result.data1.i.exp3.ICC.low[[i]],tot)
  }
}


result.data1.exp3.ICC.low<-matrix(ncol=dim(result.data1.i.exp3.ICC.low[[1]])[2])
for (i in 1:nsim){
  result.data1.exp3.ICC.low<-rbind(result.data1.exp3.ICC.low,result.data1.i.exp3.ICC.low[[i]])
}

result.data1.exp3.ICC.low<-as.data.frame(result.data1.exp3.ICC.low)
colnames(result.data1.exp3.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.ICC.low$nselected<-as.numeric(as.character(result.data1.exp3.ICC.low$nselected))
result.data1.exp3.ICC.low$tn<-as.numeric(as.character(result.data1.exp3.ICC.low$tn))
result.data1.exp3.ICC.low$fn<-as.numeric(as.character(result.data1.exp3.ICC.low$fn))
result.data1.exp3.ICC.low$fp<-as.numeric(as.character(result.data1.exp3.ICC.low$fp))
result.data1.exp3.ICC.low$tp<-as.numeric(as.character(result.data1.exp3.ICC.low$tp))

result.data1.exp3.ICC.low$sensitivity<-round(result.data1.exp3.ICC.low$tp/(result.data1.exp3.ICC.low$tp+result.data1.exp3.ICC.low$fn)*100,1)
result.data1.exp3.ICC.low$specificity<-round(result.data1.exp3.ICC.low$tn/(result.data1.exp3.ICC.low$tn+result.data1.exp3.ICC.low$fp)*100,1)
result.data1.exp3.ICC.low$false.pos.rate<-round(result.data1.exp3.ICC.low$fp/(result.data1.exp3.ICC.low$fp+result.data1.exp3.ICC.low$tn)*100,1)
result.data1.exp3.ICC.low$false.neg.rate<-round(result.data1.exp3.ICC.low$fn/(result.data1.exp3.ICC.low$fn+result.data1.exp3.ICC.low$tp)*100,1)
result.data1.exp3.ICC.low$false.disc.rate<-round(result.data1.exp3.ICC.low$fp/(result.data1.exp3.ICC.low$fp+result.data1.exp3.ICC.low$tp)*100,1)
result.data1.exp3.ICC.low$false.omit.rate<-round(result.data1.exp3.ICC.low$fn/(result.data1.exp3.ICC.low$fn+result.data1.exp3.ICC.low$tn)*100,1)


result.data1.exp3.ICC.low<-result.data1.exp3.ICC.low[-1,]
result.data1.exp3.ICC.low$sensitivity[result.data1.exp3.ICC.low$tp==0 & result.data1.exp3.ICC.low$fn>0]<-0
result.data1.exp3.ICC.low$false.disc.rate[result.data1.exp3.ICC.low$nselected==0]<-0

result.data1.exp3.ICC.low$expo<-3
result.data1.exp3.ICC.low$ICC<-0.1



# Medium ICC

RES1.i.total.exp3.ICC.med<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp3.ICC.med[[i]]<-RES1.total.exp3.ICC.med[RES1.total.exp3.ICC.med$numsim == i,]
}

result.data1.i.exp3.ICC.med<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp3.ICC.med[[i]])[2]-1)){
    RES1.i.total.exp3.ICC.med[[i]][,j]<-as.factor(as.character(RES1.i.total.exp3.ICC.med[[i]][,j]))  
    RES1.i.total.exp3.ICC.med[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp3.ICC.med[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp3.ICC.med[[i]][,j],RES1.i.total.exp3.ICC.med[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp3.ICC.med[[i]][,j],RES1.i.total.exp3.ICC.med[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp3.ICC.med[[i]][,j],RES1.i.total.exp3.ICC.med[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp3.ICC.med[[i]][,j],RES1.i.total.exp3.ICC.med[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp3.ICC.med[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp3.ICC.med[[i]]<-rbind(result.data1.i.exp3.ICC.med[[i]],tot)
  }
}


result.data1.exp3.ICC.med<-matrix(ncol=dim(result.data1.i.exp3.ICC.med[[1]])[2])
for (i in 1:nsim){
  result.data1.exp3.ICC.med<-rbind(result.data1.exp3.ICC.med,result.data1.i.exp3.ICC.med[[i]])
}

result.data1.exp3.ICC.med<-as.data.frame(result.data1.exp3.ICC.med)
colnames(result.data1.exp3.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.ICC.med$nselected<-as.numeric(as.character(result.data1.exp3.ICC.med$nselected))
result.data1.exp3.ICC.med$tn<-as.numeric(as.character(result.data1.exp3.ICC.med$tn))
result.data1.exp3.ICC.med$fn<-as.numeric(as.character(result.data1.exp3.ICC.med$fn))
result.data1.exp3.ICC.med$fp<-as.numeric(as.character(result.data1.exp3.ICC.med$fp))
result.data1.exp3.ICC.med$tp<-as.numeric(as.character(result.data1.exp3.ICC.med$tp))

result.data1.exp3.ICC.med$sensitivity<-round(result.data1.exp3.ICC.med$tp/(result.data1.exp3.ICC.med$tp+result.data1.exp3.ICC.med$fn)*100,1)
result.data1.exp3.ICC.med$specificity<-round(result.data1.exp3.ICC.med$tn/(result.data1.exp3.ICC.med$tn+result.data1.exp3.ICC.med$fp)*100,1)
result.data1.exp3.ICC.med$false.pos.rate<-round(result.data1.exp3.ICC.med$fp/(result.data1.exp3.ICC.med$fp+result.data1.exp3.ICC.med$tn)*100,1)
result.data1.exp3.ICC.med$false.neg.rate<-round(result.data1.exp3.ICC.med$fn/(result.data1.exp3.ICC.med$fn+result.data1.exp3.ICC.med$tp)*100,1)
result.data1.exp3.ICC.med$false.disc.rate<-round(result.data1.exp3.ICC.med$fp/(result.data1.exp3.ICC.med$fp+result.data1.exp3.ICC.med$tp)*100,1)
result.data1.exp3.ICC.med$false.omit.rate<-round(result.data1.exp3.ICC.med$fn/(result.data1.exp3.ICC.med$fn+result.data1.exp3.ICC.med$tn)*100,1)


result.data1.exp3.ICC.med<-result.data1.exp3.ICC.med[-1,]
result.data1.exp3.ICC.med$sensitivity[result.data1.exp3.ICC.med$tp==0 & result.data1.exp3.ICC.med$fn>0]<-0
result.data1.exp3.ICC.med$false.disc.rate[result.data1.exp3.ICC.med$nselected==0]<-0

result.data1.exp3.ICC.med$expo<-3
result.data1.exp3.ICC.med$ICC<-0.5


# High ICC
RES1.i.total.exp3.ICC.high<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp3.ICC.high[[i]]<-RES1.total.exp3.ICC.high[RES1.total.exp3.ICC.high$numsim == i,]
}

result.data1.i.exp3.ICC.high<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp3.ICC.high[[i]])[2]-1)){
    RES1.i.total.exp3.ICC.high[[i]][,j]<-as.factor(as.character(RES1.i.total.exp3.ICC.high[[i]][,j]))  
    RES1.i.total.exp3.ICC.high[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp3.ICC.high[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp3.ICC.high[[i]][,j],RES1.i.total.exp3.ICC.high[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp3.ICC.high[[i]][,j],RES1.i.total.exp3.ICC.high[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp3.ICC.high[[i]][,j],RES1.i.total.exp3.ICC.high[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp3.ICC.high[[i]][,j],RES1.i.total.exp3.ICC.high[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp3.ICC.high[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp3.ICC.high[[i]]<-rbind(result.data1.i.exp3.ICC.high[[i]],tot)
  }
}


result.data1.exp3.ICC.high<-matrix(ncol=dim(result.data1.i.exp3.ICC.high[[1]])[2])
for (i in 1:nsim){
  result.data1.exp3.ICC.high<-rbind(result.data1.exp3.ICC.high,result.data1.i.exp3.ICC.high[[i]])
}

result.data1.exp3.ICC.high<-as.data.frame(result.data1.exp3.ICC.high)
colnames(result.data1.exp3.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.ICC.high$nselected<-as.numeric(as.character(result.data1.exp3.ICC.high$nselected))
result.data1.exp3.ICC.high$tn<-as.numeric(as.character(result.data1.exp3.ICC.high$tn))
result.data1.exp3.ICC.high$fn<-as.numeric(as.character(result.data1.exp3.ICC.high$fn))
result.data1.exp3.ICC.high$fp<-as.numeric(as.character(result.data1.exp3.ICC.high$fp))
result.data1.exp3.ICC.high$tp<-as.numeric(as.character(result.data1.exp3.ICC.high$tp))

result.data1.exp3.ICC.high$sensitivity<-round(result.data1.exp3.ICC.high$tp/(result.data1.exp3.ICC.high$tp+result.data1.exp3.ICC.high$fn)*100,1)
result.data1.exp3.ICC.high$specificity<-round(result.data1.exp3.ICC.high$tn/(result.data1.exp3.ICC.high$tn+result.data1.exp3.ICC.high$fp)*100,1)
result.data1.exp3.ICC.high$false.pos.rate<-round(result.data1.exp3.ICC.high$fp/(result.data1.exp3.ICC.high$fp+result.data1.exp3.ICC.high$tn)*100,1)
result.data1.exp3.ICC.high$false.neg.rate<-round(result.data1.exp3.ICC.high$fn/(result.data1.exp3.ICC.high$fn+result.data1.exp3.ICC.high$tp)*100,1)
result.data1.exp3.ICC.high$false.disc.rate<-round(result.data1.exp3.ICC.high$fp/(result.data1.exp3.ICC.high$fp+result.data1.exp3.ICC.high$tp)*100,1)
result.data1.exp3.ICC.high$false.omit.rate<-round(result.data1.exp3.ICC.high$fn/(result.data1.exp3.ICC.high$fn+result.data1.exp3.ICC.high$tn)*100,1)


result.data1.exp3.ICC.high<-result.data1.exp3.ICC.high[-1,]
result.data1.exp3.ICC.high$sensitivity[result.data1.exp3.ICC.high$tp==0 & result.data1.exp3.ICC.high$fn>0]<-0
result.data1.exp3.ICC.high$false.disc.rate[result.data1.exp3.ICC.high$nselected==0]<-0

result.data1.exp3.ICC.high$expo<-3
result.data1.exp3.ICC.high$ICC<-0.9

############### exp5

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp5")
nexp<-5


load(file="RES1.all.exp5.2023.RData")
colnames(RES1.all.exp5)
colnames(RES1.all.exp5)<-c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES1.all.exp5<-RES1.all.exp5[,c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                                "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]

load(file="twostep/RES1avRed.all.exp5.0621.RData")
colnames(RES1avRed.all.exp5.0621)
colnames(RES1avRed.all.exp5.0621)<-c("var","numsim","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES1.DLNM.all.exp5.RData")
colnames(RES1.DLNM.all.exp5)


load(file="dlnm/RES1.DLNM.AVG.all.exp5.RData")
colnames(RES1.DLNM.AVG.all.exp5)
colnames(RES1.DLNM.AVG.all.exp5)<-c("var","numsim","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                    "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES1.total.exp5<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                        list(RES1.all.exp5,RES1avRed.all.exp5.0621,RES1.DLNM.all.exp5,RES1.DLNM.AVG.all.exp5))


### Merge avec ICC

ICC.data1.exp5<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data1.exp5)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data1.exp5<-rbind(ICC.data1.exp5,ICC.i[[i]])
}
ICC.data1.exp5<-ICC.data1.exp5[-1,]
ICC.data1.exp5$numsim<-as.integer(ICC.data1.exp5$numsim)
colnames(ICC.data1.exp5)<-c("var","numsim","ICC")

RES1.total.exp5.ICC<-merge(RES1.total.exp5,ICC.data1.exp5,by=c("var","numsim"))

table(RES1.total.exp5.ICC$true.pred,RES1.total.exp5.ICC$ICC)



RES1.total.exp5.ICC.low<-RES1.total.exp5.ICC[RES1.total.exp5.ICC$ICC==0.1,]
RES1.total.exp5.ICC.med<-RES1.total.exp5.ICC[RES1.total.exp5.ICC$ICC==0.5,]
RES1.total.exp5.ICC.high<-RES1.total.exp5.ICC[RES1.total.exp5.ICC$ICC==0.9,]


### All


RES1.i.total.exp5.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp5.ICC[[i]]<-RES1.total.exp5.ICC[RES1.total.exp5.ICC$numsim == i,]
}

result.data1.i.exp5<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp5.ICC[[i]])[2]-1)){
    RES1.i.total.exp5.ICC[[i]][,j]<-as.factor(as.character(RES1.i.total.exp5.ICC[[i]][,j]))  
    RES1.i.total.exp5.ICC[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp5.ICC[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp5.ICC[[i]][,j],RES1.i.total.exp5.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp5.ICC[[i]][,j],RES1.i.total.exp5.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp5.ICC[[i]][,j],RES1.i.total.exp5.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp5.ICC[[i]][,j],RES1.i.total.exp5.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp5.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp5[[i]]<-rbind(result.data1.i.exp5[[i]],tot)
  }
}


result.data1.exp5<-matrix(ncol=dim(result.data1.i.exp5[[1]])[2])
for (i in 1:nsim){
  result.data1.exp5<-rbind(result.data1.exp5,result.data1.i.exp5[[i]])
}


result.data1.exp5<-as.data.frame(result.data1.exp5)
colnames(result.data1.exp5)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp5$nselected<-as.numeric(as.character(result.data1.exp5$nselected))
result.data1.exp5$tn<-as.numeric(as.character(result.data1.exp5$tn))
result.data1.exp5$fn<-as.numeric(as.character(result.data1.exp5$fn))
result.data1.exp5$fp<-as.numeric(as.character(result.data1.exp5$fp))
result.data1.exp5$tp<-as.numeric(as.character(result.data1.exp5$tp))

result.data1.exp5$sensitivity<-round(result.data1.exp5$tp/(result.data1.exp5$tp+result.data1.exp5$fn)*100,1)
result.data1.exp5$specificity<-round(result.data1.exp5$tn/(result.data1.exp5$tn+result.data1.exp5$fp)*100,1)
result.data1.exp5$false.pos.rate<-round(result.data1.exp5$fp/(result.data1.exp5$fp+result.data1.exp5$tn)*100,1)
result.data1.exp5$false.neg.rate<-round(result.data1.exp5$fn/(result.data1.exp5$fn+result.data1.exp5$tp)*100,1)
result.data1.exp5$false.disc.rate<-round(result.data1.exp5$fp/(result.data1.exp5$fp+result.data1.exp5$tp)*100,1)
result.data1.exp5$false.omit.rate<-round(result.data1.exp5$fn/(result.data1.exp5$fn+result.data1.exp5$tn)*100,1)

result.data1.exp5<-result.data1.exp5[-1,]

result.data1.exp5$expo<-5
result.data1.exp5$ICC<-1



# Low ICC


RES1.i.total.exp5.ICC.low<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp5.ICC.low[[i]]<-RES1.total.exp5.ICC.low[RES1.total.exp5.ICC.low$numsim == i,]
}

result.data1.i.exp5.ICC.low<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp5.ICC.low[[i]])[2]-1)){
    RES1.i.total.exp5.ICC.low[[i]][,j]<-as.factor(as.character(RES1.i.total.exp5.ICC.low[[i]][,j]))  
    RES1.i.total.exp5.ICC.low[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp5.ICC.low[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp5.ICC.low[[i]][,j],RES1.i.total.exp5.ICC.low[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp5.ICC.low[[i]][,j],RES1.i.total.exp5.ICC.low[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp5.ICC.low[[i]][,j],RES1.i.total.exp5.ICC.low[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp5.ICC.low[[i]][,j],RES1.i.total.exp5.ICC.low[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp5.ICC.low[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp5.ICC.low[[i]]<-rbind(result.data1.i.exp5.ICC.low[[i]],tot)
  }
}


result.data1.exp5.ICC.low<-matrix(ncol=dim(result.data1.i.exp5.ICC.low[[1]])[2])
for (i in 1:nsim){
  result.data1.exp5.ICC.low<-rbind(result.data1.exp5.ICC.low,result.data1.i.exp5.ICC.low[[i]])
}

result.data1.exp5.ICC.low<-as.data.frame(result.data1.exp5.ICC.low)
colnames(result.data1.exp5.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp5.ICC.low$nselected<-as.numeric(as.character(result.data1.exp5.ICC.low$nselected))
result.data1.exp5.ICC.low$tn<-as.numeric(as.character(result.data1.exp5.ICC.low$tn))
result.data1.exp5.ICC.low$fn<-as.numeric(as.character(result.data1.exp5.ICC.low$fn))
result.data1.exp5.ICC.low$fp<-as.numeric(as.character(result.data1.exp5.ICC.low$fp))
result.data1.exp5.ICC.low$tp<-as.numeric(as.character(result.data1.exp5.ICC.low$tp))

result.data1.exp5.ICC.low$sensitivity<-round(result.data1.exp5.ICC.low$tp/(result.data1.exp5.ICC.low$tp+result.data1.exp5.ICC.low$fn)*100,1)
result.data1.exp5.ICC.low$specificity<-round(result.data1.exp5.ICC.low$tn/(result.data1.exp5.ICC.low$tn+result.data1.exp5.ICC.low$fp)*100,1)
result.data1.exp5.ICC.low$false.pos.rate<-round(result.data1.exp5.ICC.low$fp/(result.data1.exp5.ICC.low$fp+result.data1.exp5.ICC.low$tn)*100,1)
result.data1.exp5.ICC.low$false.neg.rate<-round(result.data1.exp5.ICC.low$fn/(result.data1.exp5.ICC.low$fn+result.data1.exp5.ICC.low$tp)*100,1)
result.data1.exp5.ICC.low$false.disc.rate<-round(result.data1.exp5.ICC.low$fp/(result.data1.exp5.ICC.low$fp+result.data1.exp5.ICC.low$tp)*100,1)
result.data1.exp5.ICC.low$false.omit.rate<-round(result.data1.exp5.ICC.low$fn/(result.data1.exp5.ICC.low$fn+result.data1.exp5.ICC.low$tn)*100,1)

result.data1.exp5.ICC.low<-result.data1.exp5.ICC.low[-1,]

result.data1.exp5.ICC.low$expo<-5
result.data1.exp5.ICC.low$ICC<-0.1



# Medium ICC


RES1.i.total.exp5.ICC.med<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp5.ICC.med[[i]]<-RES1.total.exp5.ICC.med[RES1.total.exp5.ICC.med$numsim == i,]
}

result.data1.i.exp5.ICC.med<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp5.ICC.med[[i]])[2]-1)){
    RES1.i.total.exp5.ICC.med[[i]][,j]<-as.factor(as.character(RES1.i.total.exp5.ICC.med[[i]][,j]))  
    RES1.i.total.exp5.ICC.med[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp5.ICC.med[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp5.ICC.med[[i]][,j],RES1.i.total.exp5.ICC.med[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp5.ICC.med[[i]][,j],RES1.i.total.exp5.ICC.med[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp5.ICC.med[[i]][,j],RES1.i.total.exp5.ICC.med[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp5.ICC.med[[i]][,j],RES1.i.total.exp5.ICC.med[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp5.ICC.med[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp5.ICC.med[[i]]<-rbind(result.data1.i.exp5.ICC.med[[i]],tot)
  }
}


result.data1.exp5.ICC.med<-matrix(ncol=dim(result.data1.i.exp5.ICC.med[[1]])[2])
for (i in 1:nsim){
  result.data1.exp5.ICC.med<-rbind(result.data1.exp5.ICC.med,result.data1.i.exp5.ICC.med[[i]])
}

result.data1.exp5.ICC.med<-as.data.frame(result.data1.exp5.ICC.med)
colnames(result.data1.exp5.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")


result.data1.exp5.ICC.med$nselected<-as.numeric(as.character(result.data1.exp5.ICC.med$nselected))
result.data1.exp5.ICC.med$tn<-as.numeric(as.character(result.data1.exp5.ICC.med$tn))
result.data1.exp5.ICC.med$fn<-as.numeric(as.character(result.data1.exp5.ICC.med$fn))
result.data1.exp5.ICC.med$fp<-as.numeric(as.character(result.data1.exp5.ICC.med$fp))
result.data1.exp5.ICC.med$tp<-as.numeric(as.character(result.data1.exp5.ICC.med$tp))

result.data1.exp5.ICC.med$sensitivity<-round(result.data1.exp5.ICC.med$tp/(result.data1.exp5.ICC.med$tp+result.data1.exp5.ICC.med$fn)*100,1)
result.data1.exp5.ICC.med$specificity<-round(result.data1.exp5.ICC.med$tn/(result.data1.exp5.ICC.med$tn+result.data1.exp5.ICC.med$fp)*100,1)
result.data1.exp5.ICC.med$false.pos.rate<-round(result.data1.exp5.ICC.med$fp/(result.data1.exp5.ICC.med$fp+result.data1.exp5.ICC.med$tn)*100,1)
result.data1.exp5.ICC.med$false.neg.rate<-round(result.data1.exp5.ICC.med$fn/(result.data1.exp5.ICC.med$fn+result.data1.exp5.ICC.med$tp)*100,1)
result.data1.exp5.ICC.med$false.disc.rate<-round(result.data1.exp5.ICC.med$fp/(result.data1.exp5.ICC.med$fp+result.data1.exp5.ICC.med$tp)*100,1)
result.data1.exp5.ICC.med$false.omit.rate<-round(result.data1.exp5.ICC.med$fn/(result.data1.exp5.ICC.med$fn+result.data1.exp5.ICC.med$tn)*100,1)

result.data1.exp5.ICC.med<-result.data1.exp5.ICC.med[-1,]

result.data1.exp5.ICC.med$expo<-5
result.data1.exp5.ICC.med$ICC<-0.5


# High ICC

RES1.i.total.exp5.ICC.high<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp5.ICC.high[[i]]<-RES1.total.exp5.ICC.high[RES1.total.exp5.ICC.high$numsim == i,]
}

result.data1.i.exp5.ICC.high<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp5.ICC.high[[i]])[2]-1)){
    RES1.i.total.exp5.ICC.high[[i]][,j]<-as.factor(as.character(RES1.i.total.exp5.ICC.high[[i]][,j]))  
    RES1.i.total.exp5.ICC.high[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp5.ICC.high[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp5.ICC.high[[i]][,j],RES1.i.total.exp5.ICC.high[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp5.ICC.high[[i]][,j],RES1.i.total.exp5.ICC.high[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp5.ICC.high[[i]][,j],RES1.i.total.exp5.ICC.high[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp5.ICC.high[[i]][,j],RES1.i.total.exp5.ICC.high[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp5.ICC.high[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp5.ICC.high[[i]]<-rbind(result.data1.i.exp5.ICC.high[[i]],tot)
  }
}


result.data1.exp5.ICC.high<-matrix(ncol=dim(result.data1.i.exp5.ICC.high[[1]])[2])
for (i in 1:nsim){
  result.data1.exp5.ICC.high<-rbind(result.data1.exp5.ICC.high,result.data1.i.exp5.ICC.high[[i]])
}

result.data1.exp5.ICC.high<-as.data.frame(result.data1.exp5.ICC.high)
colnames(result.data1.exp5.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp5.ICC.high$nselected<-as.numeric(as.character(result.data1.exp5.ICC.high$nselected))
result.data1.exp5.ICC.high$tn<-as.numeric(as.character(result.data1.exp5.ICC.high$tn))
result.data1.exp5.ICC.high$fn<-as.numeric(as.character(result.data1.exp5.ICC.high$fn))
result.data1.exp5.ICC.high$fp<-as.numeric(as.character(result.data1.exp5.ICC.high$fp))
result.data1.exp5.ICC.high$tp<-as.numeric(as.character(result.data1.exp5.ICC.high$tp))

result.data1.exp5.ICC.high$sensitivity<-round(result.data1.exp5.ICC.high$tp/(result.data1.exp5.ICC.high$tp+result.data1.exp5.ICC.high$fn)*100,1)
result.data1.exp5.ICC.high$specificity<-round(result.data1.exp5.ICC.high$tn/(result.data1.exp5.ICC.high$tn+result.data1.exp5.ICC.high$fp)*100,1)
result.data1.exp5.ICC.high$false.pos.rate<-round(result.data1.exp5.ICC.high$fp/(result.data1.exp5.ICC.high$fp+result.data1.exp5.ICC.high$tn)*100,1)
result.data1.exp5.ICC.high$false.neg.rate<-round(result.data1.exp5.ICC.high$fn/(result.data1.exp5.ICC.high$fn+result.data1.exp5.ICC.high$tp)*100,1)
result.data1.exp5.ICC.high$false.disc.rate<-round(result.data1.exp5.ICC.high$fp/(result.data1.exp5.ICC.high$fp+result.data1.exp5.ICC.high$tp)*100,1)
result.data1.exp5.ICC.high$false.omit.rate<-round(result.data1.exp5.ICC.high$fn/(result.data1.exp5.ICC.high$fn+result.data1.exp5.ICC.high$tn)*100,1)

result.data1.exp5.ICC.high<-result.data1.exp5.ICC.high[-1,]

result.data1.exp5.ICC.high$expo<-5
result.data1.exp5.ICC.high$ICC<-0.9


############### exp10

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp10")
nexp<-10


load(file="RES1.all.exp10.2023.RData")
colnames(RES1.all.exp10)
colnames(RES1.all.exp10)<-c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES1.all.exp10<-RES1.all.exp10[,c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                            "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES1avRed.all.exp10.0621.RData")
colnames(RES1avRed.all.exp10.0621)
colnames(RES1avRed.all.exp10.0621)<-c("var","numsim","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES1.DLNM.all.exp10.RData")
colnames(RES1.DLNM.all.exp10)


load(file="dlnm/RES1.DLNM.AVG.all.exp10.RData")
colnames(RES1.DLNM.AVG.all.exp10)
colnames(RES1.DLNM.AVG.all.exp10)<-c("var","numsim","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                    "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES1.total.exp10<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                        list(RES1.all.exp10,RES1avRed.all.exp10.0621,RES1.DLNM.all.exp10,RES1.DLNM.AVG.all.exp10))


### Merge avec ICC

ICC.data1.exp10<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data1.exp10)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data1.exp10<-rbind(ICC.data1.exp10,ICC.i[[i]])
}
ICC.data1.exp10<-ICC.data1.exp10[-1,]
ICC.data1.exp10$numsim<-as.integer(ICC.data1.exp10$numsim)
colnames(ICC.data1.exp10)<-c("var","numsim","ICC")

RES1.total.exp10.ICC<-merge(RES1.total.exp10,ICC.data1.exp10,by=c("var","numsim"))

table(RES1.total.exp10.ICC$true.pred,RES1.total.exp10.ICC$ICC)



RES1.total.exp10.ICC.low<-RES1.total.exp10.ICC[RES1.total.exp10.ICC$ICC==0.1,]
RES1.total.exp10.ICC.med<-RES1.total.exp10.ICC[RES1.total.exp10.ICC$ICC==0.5,]
RES1.total.exp10.ICC.high<-RES1.total.exp10.ICC[RES1.total.exp10.ICC$ICC==0.9,]


### All


RES1.i.total.exp10.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp10.ICC[[i]]<-RES1.total.exp10.ICC[RES1.total.exp10.ICC$numsim == i,]
}

result.data1.i.exp10<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp10.ICC[[i]])[2]-1)){
    RES1.i.total.exp10.ICC[[i]][,j]<-as.factor(as.character(RES1.i.total.exp10.ICC[[i]][,j]))  
    RES1.i.total.exp10.ICC[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp10.ICC[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp10.ICC[[i]][,j],RES1.i.total.exp10.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp10.ICC[[i]][,j],RES1.i.total.exp10.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp10.ICC[[i]][,j],RES1.i.total.exp10.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp10.ICC[[i]][,j],RES1.i.total.exp10.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp10.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp10[[i]]<-rbind(result.data1.i.exp10[[i]],tot)
  }
}


result.data1.exp10<-matrix(ncol=dim(result.data1.i.exp10[[1]])[2])
for (i in 1:nsim){
  result.data1.exp10<-rbind(result.data1.exp10,result.data1.i.exp10[[i]])
}

result.data1.exp10<-as.data.frame(result.data1.exp10)
colnames(result.data1.exp10)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp10$nselected<-as.numeric(as.character(result.data1.exp10$nselected))
result.data1.exp10$tn<-as.numeric(as.character(result.data1.exp10$tn))
result.data1.exp10$fn<-as.numeric(as.character(result.data1.exp10$fn))
result.data1.exp10$fp<-as.numeric(as.character(result.data1.exp10$fp))
result.data1.exp10$tp<-as.numeric(as.character(result.data1.exp10$tp))

result.data1.exp10$sensitivity<-round(result.data1.exp10$tp/(result.data1.exp10$tp+result.data1.exp10$fn)*100,1)
result.data1.exp10$specificity<-round(result.data1.exp10$tn/(result.data1.exp10$tn+result.data1.exp10$fp)*100,1)
result.data1.exp10$false.pos.rate<-round(result.data1.exp10$fp/(result.data1.exp10$fp+result.data1.exp10$tn)*100,1)
result.data1.exp10$false.neg.rate<-round(result.data1.exp10$fn/(result.data1.exp10$fn+result.data1.exp10$tp)*100,1)
result.data1.exp10$false.disc.rate<-round(result.data1.exp10$fp/(result.data1.exp10$fp+result.data1.exp10$tp)*100,1)
result.data1.exp10$false.omit.rate<-round(result.data1.exp10$fn/(result.data1.exp10$fn+result.data1.exp10$tn)*100,1)

result.data1.exp10<-result.data1.exp10[-1,]

result.data1.exp10$expo<-10
result.data1.exp10$ICC<-1




# Low ICC


RES1.i.total.exp10.ICC.low<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp10.ICC.low[[i]]<-RES1.total.exp10.ICC.low[RES1.total.exp10.ICC.low$numsim == i,]
}

result.data1.i.exp10.ICC.low<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp10.ICC.low[[i]])[2]-1)){
    RES1.i.total.exp10.ICC.low[[i]][,j]<-as.factor(as.character(RES1.i.total.exp10.ICC.low[[i]][,j]))  
    RES1.i.total.exp10.ICC.low[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp10.ICC.low[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp10.ICC.low[[i]][,j],RES1.i.total.exp10.ICC.low[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp10.ICC.low[[i]][,j],RES1.i.total.exp10.ICC.low[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp10.ICC.low[[i]][,j],RES1.i.total.exp10.ICC.low[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp10.ICC.low[[i]][,j],RES1.i.total.exp10.ICC.low[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp10.ICC.low[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp10.ICC.low[[i]]<-rbind(result.data1.i.exp10.ICC.low[[i]],tot)
  }
}


result.data1.exp10.ICC.low<-matrix(ncol=dim(result.data1.i.exp10.ICC.low[[1]])[2])
for (i in 1:nsim){
  result.data1.exp10.ICC.low<-rbind(result.data1.exp10.ICC.low,result.data1.i.exp10.ICC.low[[i]])
}

result.data1.exp10.ICC.low<-as.data.frame(result.data1.exp10.ICC.low)
colnames(result.data1.exp10.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp10.ICC.low$nselected<-as.numeric(as.character(result.data1.exp10.ICC.low$nselected))
result.data1.exp10.ICC.low$tn<-as.numeric(as.character(result.data1.exp10.ICC.low$tn))
result.data1.exp10.ICC.low$fn<-as.numeric(as.character(result.data1.exp10.ICC.low$fn))
result.data1.exp10.ICC.low$fp<-as.numeric(as.character(result.data1.exp10.ICC.low$fp))
result.data1.exp10.ICC.low$tp<-as.numeric(as.character(result.data1.exp10.ICC.low$tp))

result.data1.exp10.ICC.low$sensitivity<-round(result.data1.exp10.ICC.low$tp/(result.data1.exp10.ICC.low$tp+result.data1.exp10.ICC.low$fn)*100,1)
result.data1.exp10.ICC.low$specificity<-round(result.data1.exp10.ICC.low$tn/(result.data1.exp10.ICC.low$tn+result.data1.exp10.ICC.low$fp)*100,1)
result.data1.exp10.ICC.low$false.pos.rate<-round(result.data1.exp10.ICC.low$fp/(result.data1.exp10.ICC.low$fp+result.data1.exp10.ICC.low$tn)*100,1)
result.data1.exp10.ICC.low$false.neg.rate<-round(result.data1.exp10.ICC.low$fn/(result.data1.exp10.ICC.low$fn+result.data1.exp10.ICC.low$tp)*100,1)
result.data1.exp10.ICC.low$false.disc.rate<-round(result.data1.exp10.ICC.low$fp/(result.data1.exp10.ICC.low$fp+result.data1.exp10.ICC.low$tp)*100,1)
result.data1.exp10.ICC.low$false.omit.rate<-round(result.data1.exp10.ICC.low$fn/(result.data1.exp10.ICC.low$fn+result.data1.exp10.ICC.low$tn)*100,1)

result.data1.exp10.ICC.low<-result.data1.exp10.ICC.low[-1,]

result.data1.exp10.ICC.low$expo<-10
result.data1.exp10.ICC.low$ICC<-0.1



# Medium ICC


RES1.i.total.exp10.ICC.med<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp10.ICC.med[[i]]<-RES1.total.exp10.ICC.med[RES1.total.exp10.ICC.med$numsim == i,]
}

result.data1.i.exp10.ICC.med<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp10.ICC.med[[i]])[2]-1)){
    RES1.i.total.exp10.ICC.med[[i]][,j]<-as.factor(as.character(RES1.i.total.exp10.ICC.med[[i]][,j]))  
    RES1.i.total.exp10.ICC.med[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp10.ICC.med[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp10.ICC.med[[i]][,j],RES1.i.total.exp10.ICC.med[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp10.ICC.med[[i]][,j],RES1.i.total.exp10.ICC.med[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp10.ICC.med[[i]][,j],RES1.i.total.exp10.ICC.med[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp10.ICC.med[[i]][,j],RES1.i.total.exp10.ICC.med[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp10.ICC.med[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp10.ICC.med[[i]]<-rbind(result.data1.i.exp10.ICC.med[[i]],tot)
  }
}


result.data1.exp10.ICC.med<-matrix(ncol=dim(result.data1.i.exp10.ICC.med[[1]])[2])
for (i in 1:nsim){
  result.data1.exp10.ICC.med<-rbind(result.data1.exp10.ICC.med,result.data1.i.exp10.ICC.med[[i]])
}

result.data1.exp10.ICC.med<-as.data.frame(result.data1.exp10.ICC.med)
colnames(result.data1.exp10.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")


result.data1.exp10.ICC.med$nselected<-as.numeric(as.character(result.data1.exp10.ICC.med$nselected))
result.data1.exp10.ICC.med$tn<-as.numeric(as.character(result.data1.exp10.ICC.med$tn))
result.data1.exp10.ICC.med$fn<-as.numeric(as.character(result.data1.exp10.ICC.med$fn))
result.data1.exp10.ICC.med$fp<-as.numeric(as.character(result.data1.exp10.ICC.med$fp))
result.data1.exp10.ICC.med$tp<-as.numeric(as.character(result.data1.exp10.ICC.med$tp))

result.data1.exp10.ICC.med$sensitivity<-round(result.data1.exp10.ICC.med$tp/(result.data1.exp10.ICC.med$tp+result.data1.exp10.ICC.med$fn)*100,1)
result.data1.exp10.ICC.med$specificity<-round(result.data1.exp10.ICC.med$tn/(result.data1.exp10.ICC.med$tn+result.data1.exp10.ICC.med$fp)*100,1)
result.data1.exp10.ICC.med$false.pos.rate<-round(result.data1.exp10.ICC.med$fp/(result.data1.exp10.ICC.med$fp+result.data1.exp10.ICC.med$tn)*100,1)
result.data1.exp10.ICC.med$false.neg.rate<-round(result.data1.exp10.ICC.med$fn/(result.data1.exp10.ICC.med$fn+result.data1.exp10.ICC.med$tp)*100,1)
result.data1.exp10.ICC.med$false.disc.rate<-round(result.data1.exp10.ICC.med$fp/(result.data1.exp10.ICC.med$fp+result.data1.exp10.ICC.med$tp)*100,1)
result.data1.exp10.ICC.med$false.omit.rate<-round(result.data1.exp10.ICC.med$fn/(result.data1.exp10.ICC.med$fn+result.data1.exp10.ICC.med$tn)*100,1)

result.data1.exp10.ICC.med<-result.data1.exp10.ICC.med[-1,]

result.data1.exp10.ICC.med$expo<-10
result.data1.exp10.ICC.med$ICC<-0.5


# High ICC

RES1.i.total.exp10.ICC.high<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp10.ICC.high[[i]]<-RES1.total.exp10.ICC.high[RES1.total.exp10.ICC.high$numsim == i,]
}

result.data1.i.exp10.ICC.high<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp10.ICC.high[[i]])[2]-1)){
    RES1.i.total.exp10.ICC.high[[i]][,j]<-as.factor(as.character(RES1.i.total.exp10.ICC.high[[i]][,j]))  
    RES1.i.total.exp10.ICC.high[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp10.ICC.high[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp10.ICC.high[[i]][,j],RES1.i.total.exp10.ICC.high[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp10.ICC.high[[i]][,j],RES1.i.total.exp10.ICC.high[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp10.ICC.high[[i]][,j],RES1.i.total.exp10.ICC.high[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp10.ICC.high[[i]][,j],RES1.i.total.exp10.ICC.high[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp10.ICC.high[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp10.ICC.high[[i]]<-rbind(result.data1.i.exp10.ICC.high[[i]],tot)
  }
}


result.data1.exp10.ICC.high<-matrix(ncol=dim(result.data1.i.exp10.ICC.high[[1]])[2])
for (i in 1:nsim){
  result.data1.exp10.ICC.high<-rbind(result.data1.exp10.ICC.high,result.data1.i.exp10.ICC.high[[i]])
}

result.data1.exp10.ICC.high<-as.data.frame(result.data1.exp10.ICC.high)
colnames(result.data1.exp10.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp10.ICC.high$nselected<-as.numeric(as.character(result.data1.exp10.ICC.high$nselected))
result.data1.exp10.ICC.high$tn<-as.numeric(as.character(result.data1.exp10.ICC.high$tn))
result.data1.exp10.ICC.high$fn<-as.numeric(as.character(result.data1.exp10.ICC.high$fn))
result.data1.exp10.ICC.high$fp<-as.numeric(as.character(result.data1.exp10.ICC.high$fp))
result.data1.exp10.ICC.high$tp<-as.numeric(as.character(result.data1.exp10.ICC.high$tp))

result.data1.exp10.ICC.high$sensitivity<-round(result.data1.exp10.ICC.high$tp/(result.data1.exp10.ICC.high$tp+result.data1.exp10.ICC.high$fn)*100,1)
result.data1.exp10.ICC.high$specificity<-round(result.data1.exp10.ICC.high$tn/(result.data1.exp10.ICC.high$tn+result.data1.exp10.ICC.high$fp)*100,1)
result.data1.exp10.ICC.high$false.pos.rate<-round(result.data1.exp10.ICC.high$fp/(result.data1.exp10.ICC.high$fp+result.data1.exp10.ICC.high$tn)*100,1)
result.data1.exp10.ICC.high$false.neg.rate<-round(result.data1.exp10.ICC.high$fn/(result.data1.exp10.ICC.high$fn+result.data1.exp10.ICC.high$tp)*100,1)
result.data1.exp10.ICC.high$false.disc.rate<-round(result.data1.exp10.ICC.high$fp/(result.data1.exp10.ICC.high$fp+result.data1.exp10.ICC.high$tp)*100,1)
result.data1.exp10.ICC.high$false.omit.rate<-round(result.data1.exp10.ICC.high$fn/(result.data1.exp10.ICC.high$fn+result.data1.exp10.ICC.high$tn)*100,1)

result.data1.exp10.ICC.high<-result.data1.exp10.ICC.high[-1,]

result.data1.exp10.ICC.high$expo<-10
result.data1.exp10.ICC.high$ICC<-0.9

result.data1.ICC<-rbind(result.data1.exp3,result.data1.exp3.ICC.low,result.data1.exp3.ICC.med,result.data1.exp3.ICC.high,
                        result.data1.exp5,result.data1.exp5.ICC.low,result.data1.exp5.ICC.med,result.data1.exp5.ICC.high,
                        result.data1.exp10,result.data1.exp10.ICC.low,result.data1.exp10.ICC.med,result.data1.exp10.ICC.high)


result.data1.ICC$total.candidat<-result.data1.ICC$fp+result.data1.ICC$tp+result.data1.ICC$fn+result.data1.ICC$tn
result.data1.ICC$pct_fp<-result.data1.ICC$fp/result.data1.ICC$total.candidat*100
result.data1.ICC$pct_tp<-result.data1.ICC$tp/result.data1.ICC$total.candidat*100
result.data1.ICC$pct_fn<-result.data1.ICC$fn/result.data1.ICC$total.candidat*100
result.data1.ICC$pct_tn<-result.data1.ICC$tn/result.data1.ICC$total.candidat*100


save(result.data1.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data1_ICC_2023.Rdata")
write.csv2(result.data1.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data1_ICC_2023.csv")



data1.ICC<-result.data1.ICC %>% 
  group_by(Method,expo,ICC) %>% 
  summarise(
    pct_fp_mean = mean(pct_fp),
    pct_fp_sd = sd(pct_fp),
    pct_fn_mean = mean(pct_fn),
    pct_fn_sd = sd(pct_fn),
  )

data1.ICC<-data1.ICC[data1.ICC$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.bh","Raw.Exwas.bon","Raw.ExWAS.by"),]

ggplot(data = data1.ICC)+
  geom_point(aes(x=pct_fp_mean,y=ICC,shape=as.character(ICC)),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,25)+
  ylim(0,1)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - Raw data \nPerformance to identify the true exposure at the true time point by ICC")


test<-result.data1.ICC[result.data1.ICC$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.bh","Raw.Exwas.bon","Raw.ExWAS.by"),]

ggplot(data = result.data1.ICC, aes(x=ICC,y=pct_fp, group=as.factor(as.character(ICC))))+
  geom_boxplot()+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlab("ICC")+
  ylab("False positive rate")+
  ggtitle("Scenario 1: All time points associated with Y - Raw data \nFalse positive rate by ICC")



############################################# data2




############### EXP3

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY2andX/exp3")
nexp<-3


load(file="RES2.all.exp3.2023.RData")
colnames(RES2.all.exp3)
colnames(RES2.all.exp3)<-c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES2.all.exp3<-RES2.all.exp3[,c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES2avRed.all.exp3.0621.RData")
colnames(RES2avRed.all.exp3.0621)
colnames(RES2avRed.all.exp3.0621)<-c("var","numsim","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES2.DLNM.all.exp3.RData")
colnames(RES2.DLNM.all.exp3)


load(file="dlnm/RES2.DLNM.AVG.all.exp3.RData")
colnames(RES2.DLNM.AVG.all.exp3)
colnames(RES2.DLNM.AVG.all.exp3)<-c("var","numsim","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                    "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES2.total.exp3<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                        list(RES2.all.exp3,RES2avRed.all.exp3.0621,RES2.DLNM.all.exp3,RES2.DLNM.AVG.all.exp3))


### Merge avec ICC

ICC.data2.exp3<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data2.exp3)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data2.exp3<-rbind(ICC.data2.exp3,ICC.i[[i]])
}
ICC.data2.exp3<-ICC.data2.exp3[-1,]
ICC.data2.exp3$numsim<-as.integer(ICC.data2.exp3$numsim)
colnames(ICC.data2.exp3)<-c("var","numsim","ICC")

RES2.total.exp3.ICC<-merge(RES2.total.exp3,ICC.data2.exp3,by=c("var","numsim"))

table(RES2.total.exp3.ICC$true.pred,RES2.total.exp3.ICC$ICC)



########################## Performance to identify the true exposure at the true time point ####################################


RES2.total.exp3.ICC.low<-RES2.total.exp3.ICC[RES2.total.exp3.ICC$ICC==0.1,]
RES2.total.exp3.ICC.med<-RES2.total.exp3.ICC[RES2.total.exp3.ICC$ICC==0.5,]
RES2.total.exp3.ICC.high<-RES2.total.exp3.ICC[RES2.total.exp3.ICC$ICC==0.9,]


### All


RES2.i.total.exp3.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES2.i.total.exp3.ICC[[i]]<-RES2.total.exp3.ICC[RES2.total.exp3.ICC$numsim == i,]
}

result.data2.i.exp3<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES2.i.total.exp3.ICC[[i]])[2]-1)){
    RES2.i.total.exp3.ICC[[i]][,j]<-as.factor(as.character(RES2.i.total.exp3.ICC[[i]][,j]))  
    RES2.i.total.exp3.ICC[[i]]$true.pred<-as.factor(as.character(RES2.i.total.exp3.ICC[[i]]$true.pred))
    tn<-squareTable(RES2.i.total.exp3.ICC[[i]][,j],RES2.i.total.exp3.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES2.i.total.exp3.ICC[[i]][,j],RES2.i.total.exp3.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES2.i.total.exp3.ICC[[i]][,j],RES2.i.total.exp3.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES2.i.total.exp3.ICC[[i]][,j],RES2.i.total.exp3.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES2.i.total.exp3.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data2.i.exp3[[i]]<-rbind(result.data2.i.exp3[[i]],tot)
  }
}


result.data2.exp3<-matrix(ncol=dim(result.data2.i.exp3[[1]])[2])
for (i in 1:nsim){
  result.data2.exp3<-rbind(result.data2.exp3,result.data2.i.exp3[[i]])
}


result.data2.exp3<-as.data.frame(result.data2.exp3)
colnames(result.data2.exp3)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp3$nselected<-as.numeric(as.character(result.data2.exp3$nselected))
result.data2.exp3$tn<-as.numeric(as.character(result.data2.exp3$tn))
result.data2.exp3$fn<-as.numeric(as.character(result.data2.exp3$fn))
result.data2.exp3$fp<-as.numeric(as.character(result.data2.exp3$fp))
result.data2.exp3$tp<-as.numeric(as.character(result.data2.exp3$tp))

result.data2.exp3$sensitivity<-round(result.data2.exp3$tp/(result.data2.exp3$tp+result.data2.exp3$fn)*100,1)
result.data2.exp3$specificity<-round(result.data2.exp3$tn/(result.data2.exp3$tn+result.data2.exp3$fp)*100,1)
result.data2.exp3$false.pos.rate<-round(result.data2.exp3$fp/(result.data2.exp3$fp+result.data2.exp3$tn)*100,1)
result.data2.exp3$false.neg.rate<-round(result.data2.exp3$fn/(result.data2.exp3$fn+result.data2.exp3$tp)*100,1)
result.data2.exp3$false.disc.rate<-round(result.data2.exp3$fp/(result.data2.exp3$fp+result.data2.exp3$tp)*100,1)
result.data2.exp3$false.omit.rate<-round(result.data2.exp3$fn/(result.data2.exp3$fn+result.data2.exp3$tn)*100,1)

result.data2.exp3$expo<-3
result.data2.exp3$ICC<-1


# Low ICC
result.data2.exp3.ICC.low<-NULL
for (j in 4:(dim(RES2.total.exp3.ICC.low)[2]-1)){
  RES2.total.exp3.ICC.low[,j]<-as.factor(as.character(RES2.total.exp3.ICC.low[,j]))  
  RES2.total.exp3.ICC.low$true.pred<-as.factor(as.character(RES2.total.exp3.ICC.low$true.pred))
  tn<-squareTable(RES2.total.exp3.ICC.low[,j],RES2.total.exp3.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp3.ICC.low[,j],RES2.total.exp3.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp3.ICC.low[,j],RES2.total.exp3.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp3.ICC.low[,j],RES2.total.exp3.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp3.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp3.ICC.low<-rbind(result.data2.exp3.ICC.low,tot)
}


result.data2.exp3.ICC.low<-as.data.frame(result.data2.exp3.ICC.low)
colnames(result.data2.exp3.ICC.low)<-c("Method","nselected","tn","fn","fp","tp",'data.i')

result.data2.exp3.ICC.low$nselected<-as.numeric(as.character(result.data2.exp3.ICC.low$nselected))
result.data2.exp3.ICC.low$tn<-as.numeric(as.character(result.data2.exp3.ICC.low$tn))
result.data2.exp3.ICC.low$fn<-as.numeric(as.character(result.data2.exp3.ICC.low$fn))
result.data2.exp3.ICC.low$fp<-as.numeric(as.character(result.data2.exp3.ICC.low$fp))
result.data2.exp3.ICC.low$tp<-as.numeric(as.character(result.data2.exp3.ICC.low$tp))

result.data2.exp3.ICC.low$sensitivity<-round(result.data2.exp3.ICC.low$tp/(result.data2.exp3.ICC.low$tp+result.data2.exp3.ICC.low$fn)*100,1)
result.data2.exp3.ICC.low$specificity<-round(result.data2.exp3.ICC.low$tn/(result.data2.exp3.ICC.low$tn+result.data2.exp3.ICC.low$fp)*100,1)
result.data2.exp3.ICC.low$false.pos.rate<-round(result.data2.exp3.ICC.low$fp/(result.data2.exp3.ICC.low$fp+result.data2.exp3.ICC.low$tn)*100,1)
result.data2.exp3.ICC.low$false.neg.rate<-round(result.data2.exp3.ICC.low$fn/(result.data2.exp3.ICC.low$fn+result.data2.exp3.ICC.low$tp)*100,1)
result.data2.exp3.ICC.low$false.disc.rate<-round(result.data2.exp3.ICC.low$fp/(result.data2.exp3.ICC.low$fp+result.data2.exp3.ICC.low$tp)*100,1)
result.data2.exp3.ICC.low$false.omit.rate<-round(result.data2.exp3.ICC.low$fn/(result.data2.exp3.ICC.low$fn+result.data2.exp3.ICC.low$tn)*100,1)

result.data2.exp3.ICC.low$expo<-3
result.data2.exp3.ICC.low$ICC<-0.1



# Medium ICC
result.data2.exp3.ICC.med<-NULL
for (j in 4:(dim(RES2.total.exp3.ICC.med)[2]-1)){
  RES2.total.exp3.ICC.med[,j]<-as.factor(as.character(RES2.total.exp3.ICC.med[,j]))  
  RES2.total.exp3.ICC.med$true.pred<-as.factor(as.character(RES2.total.exp3.ICC.med$true.pred))
  tn<-squareTable(RES2.total.exp3.ICC.med[,j],RES2.total.exp3.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp3.ICC.med[,j],RES2.total.exp3.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp3.ICC.med[,j],RES2.total.exp3.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp3.ICC.med[,j],RES2.total.exp3.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp3.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp3.ICC.med<-rbind(result.data2.exp3.ICC.med,tot)
}


result.data2.exp3.ICC.med<-as.data.frame(result.data2.exp3.ICC.med)
colnames(result.data2.exp3.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp3.ICC.med$nselected<-as.numeric(as.character(result.data2.exp3.ICC.med$nselected))
result.data2.exp3.ICC.med$tn<-as.numeric(as.character(result.data2.exp3.ICC.med$tn))
result.data2.exp3.ICC.med$fn<-as.numeric(as.character(result.data2.exp3.ICC.med$fn))
result.data2.exp3.ICC.med$fp<-as.numeric(as.character(result.data2.exp3.ICC.med$fp))
result.data2.exp3.ICC.med$tp<-as.numeric(as.character(result.data2.exp3.ICC.med$tp))

result.data2.exp3.ICC.med$sensitivity<-round(result.data2.exp3.ICC.med$tp/(result.data2.exp3.ICC.med$tp+result.data2.exp3.ICC.med$fn)*100,1)
result.data2.exp3.ICC.med$specificity<-round(result.data2.exp3.ICC.med$tn/(result.data2.exp3.ICC.med$tn+result.data2.exp3.ICC.med$fp)*100,1)
result.data2.exp3.ICC.med$false.pos.rate<-round(result.data2.exp3.ICC.med$fp/(result.data2.exp3.ICC.med$fp+result.data2.exp3.ICC.med$tn)*100,1)
result.data2.exp3.ICC.med$false.neg.rate<-round(result.data2.exp3.ICC.med$fn/(result.data2.exp3.ICC.med$fn+result.data2.exp3.ICC.med$tp)*100,1)
result.data2.exp3.ICC.med$false.disc.rate<-round(result.data2.exp3.ICC.med$fp/(result.data2.exp3.ICC.med$fp+result.data2.exp3.ICC.med$tp)*100,1)
result.data2.exp3.ICC.med$false.omit.rate<-round(result.data2.exp3.ICC.med$fn/(result.data2.exp3.ICC.med$fn+result.data2.exp3.ICC.med$tn)*100,1)

result.data2.exp3.ICC.med$expo<-3
result.data2.exp3.ICC.med$ICC<-0.5


# High ICC
result.data2.exp3.ICC.high<-NULL
for (j in 4:(dim(RES2.total.exp3.ICC.high)[2]-1)){
  RES2.total.exp3.ICC.high[,j]<-as.factor(as.character(RES2.total.exp3.ICC.high[,j]))  
  RES2.total.exp3.ICC.high$true.pred<-as.factor(as.character(RES2.total.exp3.ICC.high$true.pred))
  tn<-squareTable(RES2.total.exp3.ICC.high[,j],RES2.total.exp3.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp3.ICC.high[,j],RES2.total.exp3.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp3.ICC.high[,j],RES2.total.exp3.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp3.ICC.high[,j],RES2.total.exp3.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp3.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp3.ICC.high<-rbind(result.data2.exp3.ICC.high,tot)
}


result.data2.exp3.ICC.high<-as.data.frame(result.data2.exp3.ICC.high)
colnames(result.data2.exp3.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp3.ICC.high$nselected<-as.numeric(as.character(result.data2.exp3.ICC.high$nselected))
result.data2.exp3.ICC.high$tn<-as.numeric(as.character(result.data2.exp3.ICC.high$tn))
result.data2.exp3.ICC.high$fn<-as.numeric(as.character(result.data2.exp3.ICC.high$fn))
result.data2.exp3.ICC.high$fp<-as.numeric(as.character(result.data2.exp3.ICC.high$fp))
result.data2.exp3.ICC.high$tp<-as.numeric(as.character(result.data2.exp3.ICC.high$tp))

result.data2.exp3.ICC.high$sensitivity<-round(result.data2.exp3.ICC.high$tp/(result.data2.exp3.ICC.high$tp+result.data2.exp3.ICC.high$fn)*100,1)
result.data2.exp3.ICC.high$specificity<-round(result.data2.exp3.ICC.high$tn/(result.data2.exp3.ICC.high$tn+result.data2.exp3.ICC.high$fp)*100,1)
result.data2.exp3.ICC.high$false.pos.rate<-round(result.data2.exp3.ICC.high$fp/(result.data2.exp3.ICC.high$fp+result.data2.exp3.ICC.high$tn)*100,1)
result.data2.exp3.ICC.high$false.neg.rate<-round(result.data2.exp3.ICC.high$fn/(result.data2.exp3.ICC.high$fn+result.data2.exp3.ICC.high$tp)*100,1)
result.data2.exp3.ICC.high$false.disc.rate<-round(result.data2.exp3.ICC.high$fp/(result.data2.exp3.ICC.high$fp+result.data2.exp3.ICC.high$tp)*100,1)
result.data2.exp3.ICC.high$false.omit.rate<-round(result.data2.exp3.ICC.high$fn/(result.data2.exp3.ICC.high$fn+result.data2.exp3.ICC.high$tn)*100,1)

result.data2.exp3.ICC.high$expo<-3
result.data2.exp3.ICC.high$ICC<-0.9


############### exp5

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY2andX/exp5")
nexp<-5


load(file="RES2.all.exp5.2023.RData")
colnames(RES2.all.exp5)
colnames(RES2.all.exp5)<-c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")

RES2.all.exp5<-RES2.all.exp5[,c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]

load(file="twostep/RES2avRed.all.exp5.0621.RData")
colnames(RES2avRed.all.exp5.0621)
colnames(RES2avRed.all.exp5.0621)<-c("var","numsim","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES2.DLNM.all.exp5.RData")
colnames(RES2.DLNM.all.exp5)


load(file="dlnm/RES2.DLNM.AVG.all.exp5.RData")
colnames(RES2.DLNM.AVG.all.exp5)
colnames(RES2.DLNM.AVG.all.exp5)<-c("var","numsim","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                    "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES2.total.exp5<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                        list(RES2.all.exp5,RES2avRed.all.exp5.0621,RES2.DLNM.all.exp5,RES2.DLNM.AVG.all.exp5))


### Merge avec ICC

ICC.data2.exp5<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data2.exp5)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data2.exp5<-rbind(ICC.data2.exp5,ICC.i[[i]])
}
ICC.data2.exp5<-ICC.data2.exp5[-1,]
ICC.data2.exp5$numsim<-as.integer(ICC.data2.exp5$numsim)
colnames(ICC.data2.exp5)<-c("var","numsim","ICC")

RES2.total.exp5.ICC<-merge(RES2.total.exp5,ICC.data2.exp5,by=c("var","numsim"))

table(RES2.total.exp5.ICC$true.pred,RES2.total.exp5.ICC$ICC)



RES2.total.exp5.ICC.low<-RES2.total.exp5.ICC[RES2.total.exp5.ICC$ICC==0.1,]
RES2.total.exp5.ICC.med<-RES2.total.exp5.ICC[RES2.total.exp5.ICC$ICC==0.5,]
RES2.total.exp5.ICC.high<-RES2.total.exp5.ICC[RES2.total.exp5.ICC$ICC==0.9,]


### All


RES2.i.total.exp5.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES2.i.total.exp5.ICC[[i]]<-RES2.total.exp5.ICC[RES2.total.exp5.ICC$numsim == i,]
}

result.data2.i.exp5<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES2.i.total.exp5.ICC[[i]])[2]-1)){
    RES2.i.total.exp5.ICC[[i]][,j]<-as.factor(as.character(RES2.i.total.exp5.ICC[[i]][,j]))  
    RES2.i.total.exp5.ICC[[i]]$true.pred<-as.factor(as.character(RES2.i.total.exp5.ICC[[i]]$true.pred))
    tn<-squareTable(RES2.i.total.exp5.ICC[[i]][,j],RES2.i.total.exp5.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES2.i.total.exp5.ICC[[i]][,j],RES2.i.total.exp5.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES2.i.total.exp5.ICC[[i]][,j],RES2.i.total.exp5.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES2.i.total.exp5.ICC[[i]][,j],RES2.i.total.exp5.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES2.i.total.exp5.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data2.i.exp5[[i]]<-rbind(result.data2.i.exp5[[i]],tot)
  }
}


result.data2.exp5<-matrix(ncol=dim(result.data2.i.exp5[[1]])[2])
for (i in 1:nsim){
  result.data2.exp5<-rbind(result.data2.exp5,result.data2.i.exp5[[i]])
}


result.data2.exp5<-as.data.frame(result.data2.exp5)
colnames(result.data2.exp5)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5$nselected<-as.numeric(as.character(result.data2.exp5$nselected))
result.data2.exp5$tn<-as.numeric(as.character(result.data2.exp5$tn))
result.data2.exp5$fn<-as.numeric(as.character(result.data2.exp5$fn))
result.data2.exp5$fp<-as.numeric(as.character(result.data2.exp5$fp))
result.data2.exp5$tp<-as.numeric(as.character(result.data2.exp5$tp))

result.data2.exp5$sensitivity<-round(result.data2.exp5$tp/(result.data2.exp5$tp+result.data2.exp5$fn)*100,1)
result.data2.exp5$specificity<-round(result.data2.exp5$tn/(result.data2.exp5$tn+result.data2.exp5$fp)*100,1)
result.data2.exp5$false.pos.rate<-round(result.data2.exp5$fp/(result.data2.exp5$fp+result.data2.exp5$tn)*100,1)
result.data2.exp5$false.neg.rate<-round(result.data2.exp5$fn/(result.data2.exp5$fn+result.data2.exp5$tp)*100,1)
result.data2.exp5$false.disc.rate<-round(result.data2.exp5$fp/(result.data2.exp5$fp+result.data2.exp5$tp)*100,1)
result.data2.exp5$false.omit.rate<-round(result.data2.exp5$fn/(result.data2.exp5$fn+result.data2.exp5$tn)*100,1)

result.data2.exp5$expo<-5
result.data2.exp5$ICC<-1


# Low ICC
result.data2.exp5.ICC.low<-NULL
for (j in 4:(dim(RES2.total.exp5.ICC.low)[2]-1)){
  RES2.total.exp5.ICC.low[,j]<-as.factor(as.character(RES2.total.exp5.ICC.low[,j]))  
  RES2.total.exp5.ICC.low$true.pred<-as.factor(as.character(RES2.total.exp5.ICC.low$true.pred))
  tn<-squareTable(RES2.total.exp5.ICC.low[,j],RES2.total.exp5.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp5.ICC.low[,j],RES2.total.exp5.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp5.ICC.low[,j],RES2.total.exp5.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp5.ICC.low[,j],RES2.total.exp5.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp5.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp5.ICC.low<-rbind(result.data2.exp5.ICC.low,tot)
}


result.data2.exp5.ICC.low<-as.data.frame(result.data2.exp5.ICC.low)
colnames(result.data2.exp5.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.ICC.low$nselected<-as.numeric(as.character(result.data2.exp5.ICC.low$nselected))
result.data2.exp5.ICC.low$tn<-as.numeric(as.character(result.data2.exp5.ICC.low$tn))
result.data2.exp5.ICC.low$fn<-as.numeric(as.character(result.data2.exp5.ICC.low$fn))
result.data2.exp5.ICC.low$fp<-as.numeric(as.character(result.data2.exp5.ICC.low$fp))
result.data2.exp5.ICC.low$tp<-as.numeric(as.character(result.data2.exp5.ICC.low$tp))

result.data2.exp5.ICC.low$sensitivity<-round(result.data2.exp5.ICC.low$tp/(result.data2.exp5.ICC.low$tp+result.data2.exp5.ICC.low$fn)*100,1)
result.data2.exp5.ICC.low$specificity<-round(result.data2.exp5.ICC.low$tn/(result.data2.exp5.ICC.low$tn+result.data2.exp5.ICC.low$fp)*100,1)
result.data2.exp5.ICC.low$false.pos.rate<-round(result.data2.exp5.ICC.low$fp/(result.data2.exp5.ICC.low$fp+result.data2.exp5.ICC.low$tn)*100,1)
result.data2.exp5.ICC.low$false.neg.rate<-round(result.data2.exp5.ICC.low$fn/(result.data2.exp5.ICC.low$fn+result.data2.exp5.ICC.low$tp)*100,1)
result.data2.exp5.ICC.low$false.disc.rate<-round(result.data2.exp5.ICC.low$fp/(result.data2.exp5.ICC.low$fp+result.data2.exp5.ICC.low$tp)*100,1)
result.data2.exp5.ICC.low$false.omit.rate<-round(result.data2.exp5.ICC.low$fn/(result.data2.exp5.ICC.low$fn+result.data2.exp5.ICC.low$tn)*100,1)

result.data2.exp5.ICC.low$expo<-5
result.data2.exp5.ICC.low$ICC<-0.1



# Medium ICC
result.data2.exp5.ICC.med<-NULL
for (j in 4:(dim(RES2.total.exp5.ICC.med)[2]-1)){
  RES2.total.exp5.ICC.med[,j]<-as.factor(as.character(RES2.total.exp5.ICC.med[,j]))  
  RES2.total.exp5.ICC.med$true.pred<-as.factor(as.character(RES2.total.exp5.ICC.med$true.pred))
  tn<-squareTable(RES2.total.exp5.ICC.med[,j],RES2.total.exp5.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp5.ICC.med[,j],RES2.total.exp5.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp5.ICC.med[,j],RES2.total.exp5.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp5.ICC.med[,j],RES2.total.exp5.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp5.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp5.ICC.med<-rbind(result.data2.exp5.ICC.med,tot)
}


result.data2.exp5.ICC.med<-as.data.frame(result.data2.exp5.ICC.med)
colnames(result.data2.exp5.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.ICC.med$nselected<-as.numeric(as.character(result.data2.exp5.ICC.med$nselected))
result.data2.exp5.ICC.med$tn<-as.numeric(as.character(result.data2.exp5.ICC.med$tn))
result.data2.exp5.ICC.med$fn<-as.numeric(as.character(result.data2.exp5.ICC.med$fn))
result.data2.exp5.ICC.med$fp<-as.numeric(as.character(result.data2.exp5.ICC.med$fp))
result.data2.exp5.ICC.med$tp<-as.numeric(as.character(result.data2.exp5.ICC.med$tp))

result.data2.exp5.ICC.med$sensitivity<-round(result.data2.exp5.ICC.med$tp/(result.data2.exp5.ICC.med$tp+result.data2.exp5.ICC.med$fn)*100,1)
result.data2.exp5.ICC.med$specificity<-round(result.data2.exp5.ICC.med$tn/(result.data2.exp5.ICC.med$tn+result.data2.exp5.ICC.med$fp)*100,1)
result.data2.exp5.ICC.med$false.pos.rate<-round(result.data2.exp5.ICC.med$fp/(result.data2.exp5.ICC.med$fp+result.data2.exp5.ICC.med$tn)*100,1)
result.data2.exp5.ICC.med$false.neg.rate<-round(result.data2.exp5.ICC.med$fn/(result.data2.exp5.ICC.med$fn+result.data2.exp5.ICC.med$tp)*100,1)
result.data2.exp5.ICC.med$false.disc.rate<-round(result.data2.exp5.ICC.med$fp/(result.data2.exp5.ICC.med$fp+result.data2.exp5.ICC.med$tp)*100,1)
result.data2.exp5.ICC.med$false.omit.rate<-round(result.data2.exp5.ICC.med$fn/(result.data2.exp5.ICC.med$fn+result.data2.exp5.ICC.med$tn)*100,1)

result.data2.exp5.ICC.med$expo<-5
result.data2.exp5.ICC.med$ICC<-0.5


# High ICC
result.data2.exp5.ICC.high<-NULL
for (j in 4:(dim(RES2.total.exp5.ICC.high)[2]-1)){
  RES2.total.exp5.ICC.high[,j]<-as.factor(as.character(RES2.total.exp5.ICC.high[,j]))  
  RES2.total.exp5.ICC.high$true.pred<-as.factor(as.character(RES2.total.exp5.ICC.high$true.pred))
  tn<-squareTable(RES2.total.exp5.ICC.high[,j],RES2.total.exp5.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp5.ICC.high[,j],RES2.total.exp5.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp5.ICC.high[,j],RES2.total.exp5.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp5.ICC.high[,j],RES2.total.exp5.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp5.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp5.ICC.high<-rbind(result.data2.exp5.ICC.high,tot)
}


result.data2.exp5.ICC.high<-as.data.frame(result.data2.exp5.ICC.high)
colnames(result.data2.exp5.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.ICC.high$nselected<-as.numeric(as.character(result.data2.exp5.ICC.high$nselected))
result.data2.exp5.ICC.high$tn<-as.numeric(as.character(result.data2.exp5.ICC.high$tn))
result.data2.exp5.ICC.high$fn<-as.numeric(as.character(result.data2.exp5.ICC.high$fn))
result.data2.exp5.ICC.high$fp<-as.numeric(as.character(result.data2.exp5.ICC.high$fp))
result.data2.exp5.ICC.high$tp<-as.numeric(as.character(result.data2.exp5.ICC.high$tp))

result.data2.exp5.ICC.high$sensitivity<-round(result.data2.exp5.ICC.high$tp/(result.data2.exp5.ICC.high$tp+result.data2.exp5.ICC.high$fn)*100,1)
result.data2.exp5.ICC.high$specificity<-round(result.data2.exp5.ICC.high$tn/(result.data2.exp5.ICC.high$tn+result.data2.exp5.ICC.high$fp)*100,1)
result.data2.exp5.ICC.high$false.pos.rate<-round(result.data2.exp5.ICC.high$fp/(result.data2.exp5.ICC.high$fp+result.data2.exp5.ICC.high$tn)*100,1)
result.data2.exp5.ICC.high$false.neg.rate<-round(result.data2.exp5.ICC.high$fn/(result.data2.exp5.ICC.high$fn+result.data2.exp5.ICC.high$tp)*100,1)
result.data2.exp5.ICC.high$false.disc.rate<-round(result.data2.exp5.ICC.high$fp/(result.data2.exp5.ICC.high$fp+result.data2.exp5.ICC.high$tp)*100,1)
result.data2.exp5.ICC.high$false.omit.rate<-round(result.data2.exp5.ICC.high$fn/(result.data2.exp5.ICC.high$fn+result.data2.exp5.ICC.high$tn)*100,1)

result.data2.exp5.ICC.high$expo<-5
result.data2.exp5.ICC.high$ICC<-0.9


############### exp10

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY2andX/exp10")
nexp<-10


load(file="RES2.all.exp10.2023.RData")
colnames(RES2.all.exp10)
colnames(RES2.all.exp10)<-c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                            "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES2.all.exp10<-RES2.all.exp10[,c("var","numsim","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                            "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES2avRed.all.exp10.0621.RData")
colnames(RES2avRed.all.exp10.0621)
colnames(RES2avRed.all.exp10.0621)<-c("var","numsim","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                      "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES2.DLNM.all.exp10.RData")
colnames(RES2.DLNM.all.exp10)


load(file="dlnm/RES2.DLNM.AVG.all.exp10.RData")
colnames(RES2.DLNM.AVG.all.exp10)
colnames(RES2.DLNM.AVG.all.exp10)<-c("var","numsim","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                     "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES2.total.exp10<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                         list(RES2.all.exp10,RES2avRed.all.exp10.0621,RES2.DLNM.all.exp10,RES2.DLNM.AVG.all.exp10))


### Merge avec ICC

ICC.data2.exp10<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data2.exp10)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data2.exp10<-rbind(ICC.data2.exp10,ICC.i[[i]])
}
ICC.data2.exp10<-ICC.data2.exp10[-1,]
ICC.data2.exp10$numsim<-as.integer(ICC.data2.exp10$numsim)
colnames(ICC.data2.exp10)<-c("var","numsim","ICC")

RES2.total.exp10.ICC<-merge(RES2.total.exp10,ICC.data2.exp10,by=c("var","numsim"))

table(RES2.total.exp10.ICC$true.pred,RES2.total.exp10.ICC$ICC)



RES2.total.exp10.ICC.low<-RES2.total.exp10.ICC[RES2.total.exp10.ICC$ICC==0.1,]
RES2.total.exp10.ICC.med<-RES2.total.exp10.ICC[RES2.total.exp10.ICC$ICC==0.5,]
RES2.total.exp10.ICC.high<-RES2.total.exp10.ICC[RES2.total.exp10.ICC$ICC==0.9,]


### All


RES2.i.total.exp10.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES2.i.total.exp10.ICC[[i]]<-RES2.total.exp10.ICC[RES2.total.exp10.ICC$numsim == i,]
}

result.data2.i.exp10<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES2.i.total.exp10.ICC[[i]])[2]-1)){
    RES2.i.total.exp10.ICC[[i]][,j]<-as.factor(as.character(RES2.i.total.exp10.ICC[[i]][,j]))  
    RES2.i.total.exp10.ICC[[i]]$true.pred<-as.factor(as.character(RES2.i.total.exp10.ICC[[i]]$true.pred))
    tn<-squareTable(RES2.i.total.exp10.ICC[[i]][,j],RES2.i.total.exp10.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES2.i.total.exp10.ICC[[i]][,j],RES2.i.total.exp10.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES2.i.total.exp10.ICC[[i]][,j],RES2.i.total.exp10.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES2.i.total.exp10.ICC[[i]][,j],RES2.i.total.exp10.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES2.i.total.exp10.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data2.i.exp10[[i]]<-rbind(result.data2.i.exp10[[i]],tot)
  }
}


result.data2.exp10<-matrix(ncol=dim(result.data2.i.exp10[[1]])[2])
for (i in 1:nsim){
  result.data2.exp10<-rbind(result.data2.exp10,result.data2.i.exp10[[i]])
}

result.data2.exp10<-as.data.frame(result.data2.exp10)
colnames(result.data2.exp10)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10$nselected<-as.numeric(as.character(result.data2.exp10$nselected))
result.data2.exp10$tn<-as.numeric(as.character(result.data2.exp10$tn))
result.data2.exp10$fn<-as.numeric(as.character(result.data2.exp10$fn))
result.data2.exp10$fp<-as.numeric(as.character(result.data2.exp10$fp))
result.data2.exp10$tp<-as.numeric(as.character(result.data2.exp10$tp))

result.data2.exp10$sensitivity<-round(result.data2.exp10$tp/(result.data2.exp10$tp+result.data2.exp10$fn)*100,1)
result.data2.exp10$specificity<-round(result.data2.exp10$tn/(result.data2.exp10$tn+result.data2.exp10$fp)*100,1)
result.data2.exp10$false.pos.rate<-round(result.data2.exp10$fp/(result.data2.exp10$fp+result.data2.exp10$tn)*100,1)
result.data2.exp10$false.neg.rate<-round(result.data2.exp10$fn/(result.data2.exp10$fn+result.data2.exp10$tp)*100,1)
result.data2.exp10$false.disc.rate<-round(result.data2.exp10$fp/(result.data2.exp10$fp+result.data2.exp10$tp)*100,1)
result.data2.exp10$false.omit.rate<-round(result.data2.exp10$fn/(result.data2.exp10$fn+result.data2.exp10$tn)*100,1)

result.data2.exp10$expo<-10
result.data2.exp10$ICC<-1

# Low ICC
result.data2.exp10.ICC.low<-NULL
for (j in 4:(dim(RES2.total.exp10.ICC.low)[2]-1)){
  RES2.total.exp10.ICC.low[,j]<-as.factor(as.character(RES2.total.exp10.ICC.low[,j]))  
  RES2.total.exp10.ICC.low$true.pred<-as.factor(as.character(RES2.total.exp10.ICC.low$true.pred))
  tn<-squareTable(RES2.total.exp10.ICC.low[,j],RES2.total.exp10.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp10.ICC.low[,j],RES2.total.exp10.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp10.ICC.low[,j],RES2.total.exp10.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp10.ICC.low[,j],RES2.total.exp10.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp10.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp10.ICC.low<-rbind(result.data2.exp10.ICC.low,tot)
}


result.data2.exp10.ICC.low<-as.data.frame(result.data2.exp10.ICC.low)
colnames(result.data2.exp10.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10.ICC.low$nselected<-as.numeric(as.character(result.data2.exp10.ICC.low$nselected))
result.data2.exp10.ICC.low$tn<-as.numeric(as.character(result.data2.exp10.ICC.low$tn))
result.data2.exp10.ICC.low$fn<-as.numeric(as.character(result.data2.exp10.ICC.low$fn))
result.data2.exp10.ICC.low$fp<-as.numeric(as.character(result.data2.exp10.ICC.low$fp))
result.data2.exp10.ICC.low$tp<-as.numeric(as.character(result.data2.exp10.ICC.low$tp))

result.data2.exp10.ICC.low$sensitivity<-round(result.data2.exp10.ICC.low$tp/(result.data2.exp10.ICC.low$tp+result.data2.exp10.ICC.low$fn)*100,1)
result.data2.exp10.ICC.low$specificity<-round(result.data2.exp10.ICC.low$tn/(result.data2.exp10.ICC.low$tn+result.data2.exp10.ICC.low$fp)*100,1)
result.data2.exp10.ICC.low$false.pos.rate<-round(result.data2.exp10.ICC.low$fp/(result.data2.exp10.ICC.low$fp+result.data2.exp10.ICC.low$tn)*100,1)
result.data2.exp10.ICC.low$false.neg.rate<-round(result.data2.exp10.ICC.low$fn/(result.data2.exp10.ICC.low$fn+result.data2.exp10.ICC.low$tp)*100,1)
result.data2.exp10.ICC.low$false.disc.rate<-round(result.data2.exp10.ICC.low$fp/(result.data2.exp10.ICC.low$fp+result.data2.exp10.ICC.low$tp)*100,1)
result.data2.exp10.ICC.low$false.omit.rate<-round(result.data2.exp10.ICC.low$fn/(result.data2.exp10.ICC.low$fn+result.data2.exp10.ICC.low$tn)*100,1)

result.data2.exp10.ICC.low$expo<-10
result.data2.exp10.ICC.low$ICC<-0.1

# Medium ICC
result.data2.exp10.ICC.med<-NULL
for (j in 4:(dim(RES2.total.exp10.ICC.med)[2]-1)){
  RES2.total.exp10.ICC.med[,j]<-as.factor(as.character(RES2.total.exp10.ICC.med[,j]))  
  RES2.total.exp10.ICC.med$true.pred<-as.factor(as.character(RES2.total.exp10.ICC.med$true.pred))
  tn<-squareTable(RES2.total.exp10.ICC.med[,j],RES2.total.exp10.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp10.ICC.med[,j],RES2.total.exp10.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp10.ICC.med[,j],RES2.total.exp10.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp10.ICC.med[,j],RES2.total.exp10.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp10.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp10.ICC.med<-rbind(result.data2.exp10.ICC.med,tot)
}


result.data2.exp10.ICC.med<-as.data.frame(result.data2.exp10.ICC.med)
colnames(result.data2.exp10.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10.ICC.med$nselected<-as.numeric(as.character(result.data2.exp10.ICC.med$nselected))
result.data2.exp10.ICC.med$tn<-as.numeric(as.character(result.data2.exp10.ICC.med$tn))
result.data2.exp10.ICC.med$fn<-as.numeric(as.character(result.data2.exp10.ICC.med$fn))
result.data2.exp10.ICC.med$fp<-as.numeric(as.character(result.data2.exp10.ICC.med$fp))
result.data2.exp10.ICC.med$tp<-as.numeric(as.character(result.data2.exp10.ICC.med$tp))

result.data2.exp10.ICC.med$sensitivity<-round(result.data2.exp10.ICC.med$tp/(result.data2.exp10.ICC.med$tp+result.data2.exp10.ICC.med$fn)*100,1)
result.data2.exp10.ICC.med$specificity<-round(result.data2.exp10.ICC.med$tn/(result.data2.exp10.ICC.med$tn+result.data2.exp10.ICC.med$fp)*100,1)
result.data2.exp10.ICC.med$false.pos.rate<-round(result.data2.exp10.ICC.med$fp/(result.data2.exp10.ICC.med$fp+result.data2.exp10.ICC.med$tn)*100,1)
result.data2.exp10.ICC.med$false.neg.rate<-round(result.data2.exp10.ICC.med$fn/(result.data2.exp10.ICC.med$fn+result.data2.exp10.ICC.med$tp)*100,1)
result.data2.exp10.ICC.med$false.disc.rate<-round(result.data2.exp10.ICC.med$fp/(result.data2.exp10.ICC.med$fp+result.data2.exp10.ICC.med$tp)*100,1)
result.data2.exp10.ICC.med$false.omit.rate<-round(result.data2.exp10.ICC.med$fn/(result.data2.exp10.ICC.med$fn+result.data2.exp10.ICC.med$tn)*100,1)

result.data2.exp10.ICC.med$expo<-10
result.data2.exp10.ICC.med$ICC<-0.5



# High ICC
result.data2.exp10.ICC.high<-NULL
for (j in 4:(dim(RES2.total.exp10.ICC.high)[2]-1)){
  RES2.total.exp10.ICC.high[,j]<-as.factor(as.character(RES2.total.exp10.ICC.high[,j]))  
  RES2.total.exp10.ICC.high$true.pred<-as.factor(as.character(RES2.total.exp10.ICC.high$true.pred))
  tn<-squareTable(RES2.total.exp10.ICC.high[,j],RES2.total.exp10.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp10.ICC.high[,j],RES2.total.exp10.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp10.ICC.high[,j],RES2.total.exp10.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp10.ICC.high[,j],RES2.total.exp10.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp10.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp10.ICC.high<-rbind(result.data2.exp10.ICC.high,tot)
}


result.data2.exp10.ICC.high<-as.data.frame(result.data2.exp10.ICC.high)
colnames(result.data2.exp10.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10.ICC.high$nselected<-as.numeric(as.character(result.data2.exp10.ICC.high$nselected))
result.data2.exp10.ICC.high$tn<-as.numeric(as.character(result.data2.exp10.ICC.high$tn))
result.data2.exp10.ICC.high$fn<-as.numeric(as.character(result.data2.exp10.ICC.high$fn))
result.data2.exp10.ICC.high$fp<-as.numeric(as.character(result.data2.exp10.ICC.high$fp))
result.data2.exp10.ICC.high$tp<-as.numeric(as.character(result.data2.exp10.ICC.high$tp))

result.data2.exp10.ICC.high$sensitivity<-round(result.data2.exp10.ICC.high$tp/(result.data2.exp10.ICC.high$tp+result.data2.exp10.ICC.high$fn)*100,1)
result.data2.exp10.ICC.high$specificity<-round(result.data2.exp10.ICC.high$tn/(result.data2.exp10.ICC.high$tn+result.data2.exp10.ICC.high$fp)*100,1)
result.data2.exp10.ICC.high$false.pos.rate<-round(result.data2.exp10.ICC.high$fp/(result.data2.exp10.ICC.high$fp+result.data2.exp10.ICC.high$tn)*100,1)
result.data2.exp10.ICC.high$false.neg.rate<-round(result.data2.exp10.ICC.high$fn/(result.data2.exp10.ICC.high$fn+result.data2.exp10.ICC.high$tp)*100,1)
result.data2.exp10.ICC.high$false.disc.rate<-round(result.data2.exp10.ICC.high$fp/(result.data2.exp10.ICC.high$fp+result.data2.exp10.ICC.high$tp)*100,1)
result.data2.exp10.ICC.high$false.omit.rate<-round(result.data2.exp10.ICC.high$fn/(result.data2.exp10.ICC.high$fn+result.data2.exp10.ICC.high$tn)*100,1)

result.data2.exp10.ICC.high$expo<-10
result.data2.exp10.ICC.high$ICC<-0.9

result.data2.ICC<-rbind(result.data2.exp3,result.data2.exp3.ICC.low,result.data2.exp3.ICC.med,result.data2.exp3.ICC.high,
                        result.data2.exp5,result.data2.exp5.ICC.low,result.data2.exp5.ICC.med,result.data2.exp5.ICC.high,
                        result.data2.exp10,result.data2.exp10.ICC.low,result.data2.exp10.ICC.med,result.data2.exp10.ICC.high)

result.data2.ICC$false.disc.rate[result.data2.ICC$nselected==0]<-0

save(result.data2.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data2_ICC.2023.Rdata")
write.csv2(result.data2.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data2_ICC.2023.csv")







########################## Performance to identify the true exposure whatever the true time point ####################################


############################################# DATA1

############### EXP3

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp3")
nexp<-3


load(file="RES1.all.exp3.100.2023.RData")
colnames(RES1.all.exp3.100)
colnames(RES1.all.exp3.100)<-c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                           "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES1.all.exp3.100<-RES1.all.exp3.100[,c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES1avRed.all.exp3.100.RData")
colnames(RES1avRed.all.exp3.100)
colnames(RES1avRed.all.exp3.100)<-c("numsim","var","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES1.DLNM.all.exp3.100.RData")
colnames(RES1.DLNM.all.exp3.100)
colnames(RES1.DLNM.all.exp3.100)<-c("numsim","var","true.pred","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                    "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                    "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",     
                                    "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")


load(file="dlnm/RES1.DLNM.AVG.all.exp3.100.RData")
colnames(RES1.DLNM.AVG.all.exp3.100)
colnames(RES1.DLNM.AVG.all.exp3.100)<-c("numsim","var","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                    "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES1.total.exp3.100<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                        list(RES1.all.exp3.100,RES1avRed.all.exp3.100,RES1.DLNM.all.exp3.100,RES1.DLNM.AVG.all.exp3.100))


### Merge avec ICC

ICC.data1.exp3.red<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data1.exp3.red)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data1.exp3.red<-rbind(ICC.data1.exp3.red,ICC.i[[i]])
}
ICC.data1.exp3.red<-ICC.data1.exp3.red[-1,]
ICC.data1.exp3.red$numsim<-as.integer(ICC.data1.exp3.red$numsim)
colnames(ICC.data1.exp3.red)<-c("var","numsim","ICC")

ICC.data1.exp3.red$expo_name<-ifelse(nchar(as.character(ICC.data1.exp3.red$var))==4,substring(ICC.data1.exp3.red$var,first=1,last=2),
                                          ifelse(nchar(as.character(ICC.data1.exp3.red$var))==5,substring(ICC.data1.exp3.red$var,first=1,last=3),
                                                 substring(ICC.data1.exp3.red$var,first=1,last=4)))

ICC.data1.exp3.red<-ICC.data1.exp3.red[,c("expo_name","numsim","ICC")]
ICC.data1.exp3.red<-unique(ICC.data1.exp3.red)
colnames(ICC.data1.exp3.red)<-c("var","numsim","ICC")
  
RES1.total.exp3.100.ICC<-merge(RES1.total.exp3.100,ICC.data1.exp3.red,by=c("var","numsim"))

table(RES1.total.exp3.100.ICC$true.pred,RES1.total.exp3.100.ICC$ICC)


squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


RES1.total.exp3.100.ICC.low<-RES1.total.exp3.100.ICC[RES1.total.exp3.100.ICC$ICC==0.1,]
RES1.total.exp3.100.ICC.med<-RES1.total.exp3.100.ICC[RES1.total.exp3.100.ICC$ICC==0.5,]
RES1.total.exp3.100.ICC.high<-RES1.total.exp3.100.ICC[RES1.total.exp3.100.ICC$ICC==0.9,]


### All

RES1.i.total.exp3.100.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp3.100.ICC[[i]]<-RES1.total.exp3.100.ICC[RES1.total.exp3.100.ICC$numsim == i,]
}

result.data1.i.exp3.100<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp3.100.ICC[[i]])[2]-1)){
    RES1.i.total.exp3.100.ICC[[i]][,j]<-as.factor(as.character(RES1.i.total.exp3.100.ICC[[i]][,j]))  
    RES1.i.total.exp3.100.ICC[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp3.100.ICC[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp3.100.ICC[[i]][,j],RES1.i.total.exp3.100.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp3.100.ICC[[i]][,j],RES1.i.total.exp3.100.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp3.100.ICC[[i]][,j],RES1.i.total.exp3.100.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp3.100.ICC[[i]][,j],RES1.i.total.exp3.100.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp3.100.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp3.100[[i]]<-rbind(result.data1.i.exp3.100[[i]],tot)
  }
}


result.data1.exp3.100<-matrix(ncol=dim(result.data1.i.exp3.100[[1]])[2])
for (i in 1:nsim){
  result.data1.exp3.100<-rbind(result.data1.exp3.100,result.data1.i.exp3.100[[i]])
}

result.data1.exp3.100<-as.data.frame(result.data1.exp3.100)
colnames(result.data1.exp3.100)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.100$nselected<-as.numeric(as.character(result.data1.exp3.100$nselected))
result.data1.exp3.100$tn<-as.numeric(as.character(result.data1.exp3.100$tn))
result.data1.exp3.100$fn<-as.numeric(as.character(result.data1.exp3.100$fn))
result.data1.exp3.100$fp<-as.numeric(as.character(result.data1.exp3.100$fp))
result.data1.exp3.100$tp<-as.numeric(as.character(result.data1.exp3.100$tp))

result.data1.exp3.100$sensitivity<-round(result.data1.exp3.100$tp/(result.data1.exp3.100$tp+result.data1.exp3.100$fn)*100,1)
result.data1.exp3.100$specificity<-round(result.data1.exp3.100$tn/(result.data1.exp3.100$tn+result.data1.exp3.100$fp)*100,1)
result.data1.exp3.100$false.pos.rate<-round(result.data1.exp3.100$fp/(result.data1.exp3.100$fp+result.data1.exp3.100$tn)*100,1)
result.data1.exp3.100$false.neg.rate<-round(result.data1.exp3.100$fn/(result.data1.exp3.100$fn+result.data1.exp3.100$tp)*100,1)
result.data1.exp3.100$false.disc.rate<-round(result.data1.exp3.100$fp/(result.data1.exp3.100$fp+result.data1.exp3.100$tp)*100,1)
result.data1.exp3.100$false.omit.rate<-round(result.data1.exp3.100$fn/(result.data1.exp3.100$fn+result.data1.exp3.100$tn)*100,1)

result.data1.exp3.100<-result.data1.exp3.100[-1,]

result.data1.exp3.100$expo<-3
result.data1.exp3.100$ICC<-1


# Low ICC
result.data1.exp3.100.ICC.low<-NULL
for (j in 4:(dim(RES1.total.exp3.100.ICC.low)[2]-1)){
  RES1.total.exp3.100.ICC.low[,j]<-as.factor(as.character(RES1.total.exp3.100.ICC.low[,j]))  
  RES1.total.exp3.100.ICC.low$true.pred<-as.factor(as.character(RES1.total.exp3.100.ICC.low$true.pred))
  tn<-squareTable(RES1.total.exp3.100.ICC.low[,j],RES1.total.exp3.100.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp3.100.ICC.low[,j],RES1.total.exp3.100.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp3.100.ICC.low[,j],RES1.total.exp3.100.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp3.100.ICC.low[,j],RES1.total.exp3.100.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp3.100.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp3.100.ICC.low<-rbind(result.data1.exp3.100.ICC.low,tot)
}


result.data1.exp3.100.ICC.low<-as.data.frame(result.data1.exp3.100.ICC.low)
colnames(result.data1.exp3.100.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.100.ICC.low$nselected<-as.numeric(as.character(result.data1.exp3.100.ICC.low$nselected))
result.data1.exp3.100.ICC.low$tn<-as.numeric(as.character(result.data1.exp3.100.ICC.low$tn))
result.data1.exp3.100.ICC.low$fn<-as.numeric(as.character(result.data1.exp3.100.ICC.low$fn))
result.data1.exp3.100.ICC.low$fp<-as.numeric(as.character(result.data1.exp3.100.ICC.low$fp))
result.data1.exp3.100.ICC.low$tp<-as.numeric(as.character(result.data1.exp3.100.ICC.low$tp))

result.data1.exp3.100.ICC.low$sensitivity<-round(result.data1.exp3.100.ICC.low$tp/(result.data1.exp3.100.ICC.low$tp+result.data1.exp3.100.ICC.low$fn)*100,1)
result.data1.exp3.100.ICC.low$specificity<-round(result.data1.exp3.100.ICC.low$tn/(result.data1.exp3.100.ICC.low$tn+result.data1.exp3.100.ICC.low$fp)*100,1)
result.data1.exp3.100.ICC.low$false.pos.rate<-round(result.data1.exp3.100.ICC.low$fp/(result.data1.exp3.100.ICC.low$fp+result.data1.exp3.100.ICC.low$tn)*100,1)
result.data1.exp3.100.ICC.low$false.neg.rate<-round(result.data1.exp3.100.ICC.low$fn/(result.data1.exp3.100.ICC.low$fn+result.data1.exp3.100.ICC.low$tp)*100,1)
result.data1.exp3.100.ICC.low$false.disc.rate<-round(result.data1.exp3.100.ICC.low$fp/(result.data1.exp3.100.ICC.low$fp+result.data1.exp3.100.ICC.low$tp)*100,1)
result.data1.exp3.100.ICC.low$false.omit.rate<-round(result.data1.exp3.100.ICC.low$fn/(result.data1.exp3.100.ICC.low$fn+result.data1.exp3.100.ICC.low$tn)*100,1)

result.data1.exp3.100.ICC.low$expo<-3
result.data1.exp3.100.ICC.low$ICC<-0.1



# Medium ICC
result.data1.exp3.100.ICC.med<-NULL
for (j in 4:(dim(RES1.total.exp3.100.ICC.med)[2]-1)){
  RES1.total.exp3.100.ICC.med[,j]<-as.factor(as.character(RES1.total.exp3.100.ICC.med[,j]))  
  RES1.total.exp3.100.ICC.med$true.pred<-as.factor(as.character(RES1.total.exp3.100.ICC.med$true.pred))
  tn<-squareTable(RES1.total.exp3.100.ICC.med[,j],RES1.total.exp3.100.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp3.100.ICC.med[,j],RES1.total.exp3.100.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp3.100.ICC.med[,j],RES1.total.exp3.100.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp3.100.ICC.med[,j],RES1.total.exp3.100.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp3.100.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp3.100.ICC.med<-rbind(result.data1.exp3.100.ICC.med,tot)
}


result.data1.exp3.100.ICC.med<-as.data.frame(result.data1.exp3.100.ICC.med)
colnames(result.data1.exp3.100.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.100.ICC.med$nselected<-as.numeric(as.character(result.data1.exp3.100.ICC.med$nselected))
result.data1.exp3.100.ICC.med$tn<-as.numeric(as.character(result.data1.exp3.100.ICC.med$tn))
result.data1.exp3.100.ICC.med$fn<-as.numeric(as.character(result.data1.exp3.100.ICC.med$fn))
result.data1.exp3.100.ICC.med$fp<-as.numeric(as.character(result.data1.exp3.100.ICC.med$fp))
result.data1.exp3.100.ICC.med$tp<-as.numeric(as.character(result.data1.exp3.100.ICC.med$tp))

result.data1.exp3.100.ICC.med$sensitivity<-round(result.data1.exp3.100.ICC.med$tp/(result.data1.exp3.100.ICC.med$tp+result.data1.exp3.100.ICC.med$fn)*100,1)
result.data1.exp3.100.ICC.med$specificity<-round(result.data1.exp3.100.ICC.med$tn/(result.data1.exp3.100.ICC.med$tn+result.data1.exp3.100.ICC.med$fp)*100,1)
result.data1.exp3.100.ICC.med$false.pos.rate<-round(result.data1.exp3.100.ICC.med$fp/(result.data1.exp3.100.ICC.med$fp+result.data1.exp3.100.ICC.med$tn)*100,1)
result.data1.exp3.100.ICC.med$false.neg.rate<-round(result.data1.exp3.100.ICC.med$fn/(result.data1.exp3.100.ICC.med$fn+result.data1.exp3.100.ICC.med$tp)*100,1)
result.data1.exp3.100.ICC.med$false.disc.rate<-round(result.data1.exp3.100.ICC.med$fp/(result.data1.exp3.100.ICC.med$fp+result.data1.exp3.100.ICC.med$tp)*100,1)
result.data1.exp3.100.ICC.med$false.omit.rate<-round(result.data1.exp3.100.ICC.med$fn/(result.data1.exp3.100.ICC.med$fn+result.data1.exp3.100.ICC.med$tn)*100,1)

result.data1.exp3.100.ICC.med$expo<-3
result.data1.exp3.100.ICC.med$ICC<-0.5


# High ICC
result.data1.exp3.100.ICC.high<-NULL
for (j in 4:(dim(RES1.total.exp3.100.ICC.high)[2]-1)){
  RES1.total.exp3.100.ICC.high[,j]<-as.factor(as.character(RES1.total.exp3.100.ICC.high[,j]))  
  RES1.total.exp3.100.ICC.high$true.pred<-as.factor(as.character(RES1.total.exp3.100.ICC.high$true.pred))
  tn<-squareTable(RES1.total.exp3.100.ICC.high[,j],RES1.total.exp3.100.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp3.100.ICC.high[,j],RES1.total.exp3.100.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp3.100.ICC.high[,j],RES1.total.exp3.100.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp3.100.ICC.high[,j],RES1.total.exp3.100.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp3.100.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp3.100.ICC.high<-rbind(result.data1.exp3.100.ICC.high,tot)
}


result.data1.exp3.100.ICC.high<-as.data.frame(result.data1.exp3.100.ICC.high)
colnames(result.data1.exp3.100.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp3.100.ICC.high$nselected<-as.numeric(as.character(result.data1.exp3.100.ICC.high$nselected))
result.data1.exp3.100.ICC.high$tn<-as.numeric(as.character(result.data1.exp3.100.ICC.high$tn))
result.data1.exp3.100.ICC.high$fn<-as.numeric(as.character(result.data1.exp3.100.ICC.high$fn))
result.data1.exp3.100.ICC.high$fp<-as.numeric(as.character(result.data1.exp3.100.ICC.high$fp))
result.data1.exp3.100.ICC.high$tp<-as.numeric(as.character(result.data1.exp3.100.ICC.high$tp))

result.data1.exp3.100.ICC.high$sensitivity<-round(result.data1.exp3.100.ICC.high$tp/(result.data1.exp3.100.ICC.high$tp+result.data1.exp3.100.ICC.high$fn)*100,1)
result.data1.exp3.100.ICC.high$specificity<-round(result.data1.exp3.100.ICC.high$tn/(result.data1.exp3.100.ICC.high$tn+result.data1.exp3.100.ICC.high$fp)*100,1)
result.data1.exp3.100.ICC.high$false.pos.rate<-round(result.data1.exp3.100.ICC.high$fp/(result.data1.exp3.100.ICC.high$fp+result.data1.exp3.100.ICC.high$tn)*100,1)
result.data1.exp3.100.ICC.high$false.neg.rate<-round(result.data1.exp3.100.ICC.high$fn/(result.data1.exp3.100.ICC.high$fn+result.data1.exp3.100.ICC.high$tp)*100,1)
result.data1.exp3.100.ICC.high$false.disc.rate<-round(result.data1.exp3.100.ICC.high$fp/(result.data1.exp3.100.ICC.high$fp+result.data1.exp3.100.ICC.high$tp)*100,1)
result.data1.exp3.100.ICC.high$false.omit.rate<-round(result.data1.exp3.100.ICC.high$fn/(result.data1.exp3.100.ICC.high$fn+result.data1.exp3.100.ICC.high$tn)*100,1)

result.data1.exp3.100.ICC.high$expo<-3
result.data1.exp3.100.ICC.high$ICC<-0.9


############### exp5

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp5")
nexp<-5


load(file="RES1.all.exp5.100.2023.RData")
colnames(RES1.all.exp5.100)
colnames(RES1.all.exp5.100)<-c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES1.all.exp5.100<-RES1.all.exp5.100[,c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES1avRed.all.exp5.100.RData")
colnames(RES1avRed.all.exp5.100)
colnames(RES1avRed.all.exp5.100)<-c("numsim","var","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                    "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES1.DLNM.all.exp5.100.RData")
colnames(RES1.DLNM.all.exp5.100)
colnames(RES1.DLNM.all.exp5.100)<-c("numsim","var","true.pred","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                    "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                    "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",     
                                    "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")


load(file="dlnm/RES1.DLNM.AVG.all.exp5.100.RData")
colnames(RES1.DLNM.AVG.all.exp5.100)
colnames(RES1.DLNM.AVG.all.exp5.100)<-c("numsim","var","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                        "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES1.total.exp5.100<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                            list(RES1.all.exp5.100,RES1avRed.all.exp5.100,RES1.DLNM.all.exp5.100,RES1.DLNM.AVG.all.exp5.100))


### Merge avec ICC

ICC.data1.exp5.red<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data1.exp5.red)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data1.exp5.red<-rbind(ICC.data1.exp5.red,ICC.i[[i]])
}
ICC.data1.exp5.red<-ICC.data1.exp5.red[-1,]
ICC.data1.exp5.red$numsim<-as.integer(ICC.data1.exp5.red$numsim)
colnames(ICC.data1.exp5.red)<-c("var","numsim","ICC")

ICC.data1.exp5.red$expo_name<-ifelse(nchar(as.character(ICC.data1.exp5.red$var))==4,substring(ICC.data1.exp5.red$var,first=1,last=2),
                                     ifelse(nchar(as.character(ICC.data1.exp5.red$var))==5,substring(ICC.data1.exp5.red$var,first=1,last=3),
                                            substring(ICC.data1.exp5.red$var,first=1,last=4)))

ICC.data1.exp5.red<-ICC.data1.exp5.red[,c("expo_name","numsim","ICC")]
ICC.data1.exp5.red<-unique(ICC.data1.exp5.red)
colnames(ICC.data1.exp5.red)<-c("var","numsim","ICC")

RES1.total.exp5.100.ICC<-merge(RES1.total.exp5.100,ICC.data1.exp5.red,by=c("var","numsim"))

table(RES1.total.exp5.100.ICC$true.pred,RES1.total.exp5.100.ICC$ICC)


squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


RES1.total.exp5.100.ICC.low<-RES1.total.exp5.100.ICC[RES1.total.exp5.100.ICC$ICC==0.1,]
RES1.total.exp5.100.ICC.med<-RES1.total.exp5.100.ICC[RES1.total.exp5.100.ICC$ICC==0.5,]
RES1.total.exp5.100.ICC.high<-RES1.total.exp5.100.ICC[RES1.total.exp5.100.ICC$ICC==0.9,]


### All
RES1.i.total.exp5.100.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp5.100.ICC[[i]]<-RES1.total.exp5.100.ICC[RES1.total.exp5.100.ICC$numsim == i,]
}

result.data1.i.exp5.100<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp5.100.ICC[[i]])[2]-1)){
    RES1.i.total.exp5.100.ICC[[i]][,j]<-as.factor(as.character(RES1.i.total.exp5.100.ICC[[i]][,j]))  
    RES1.i.total.exp5.100.ICC[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp5.100.ICC[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp5.100.ICC[[i]][,j],RES1.i.total.exp5.100.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp5.100.ICC[[i]][,j],RES1.i.total.exp5.100.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp5.100.ICC[[i]][,j],RES1.i.total.exp5.100.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp5.100.ICC[[i]][,j],RES1.i.total.exp5.100.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp5.100.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp5.100[[i]]<-rbind(result.data1.i.exp5.100[[i]],tot)
  }
}


result.data1.exp5.100<-matrix(ncol=dim(result.data1.i.exp5.100[[1]])[2])
for (i in 1:nsim){
  result.data1.exp5.100<-rbind(result.data1.exp5.100,result.data1.i.exp5.100[[i]])
}

result.data1.exp5.100<-as.data.frame(result.data1.exp5.100)
colnames(result.data1.exp5.100)<-c("Method","nselected","tn","fn","fp","tp","data.i")


result.data1.exp5.100$nselected<-as.numeric(as.character(result.data1.exp5.100$nselected))
result.data1.exp5.100$tn<-as.numeric(as.character(result.data1.exp5.100$tn))
result.data1.exp5.100$fn<-as.numeric(as.character(result.data1.exp5.100$fn))
result.data1.exp5.100$fp<-as.numeric(as.character(result.data1.exp5.100$fp))
result.data1.exp5.100$tp<-as.numeric(as.character(result.data1.exp5.100$tp))

result.data1.exp5.100$sensitivity<-round(result.data1.exp5.100$tp/(result.data1.exp5.100$tp+result.data1.exp5.100$fn)*100,1)
result.data1.exp5.100$specificity<-round(result.data1.exp5.100$tn/(result.data1.exp5.100$tn+result.data1.exp5.100$fp)*100,1)
result.data1.exp5.100$false.pos.rate<-round(result.data1.exp5.100$fp/(result.data1.exp5.100$fp+result.data1.exp5.100$tn)*100,1)
result.data1.exp5.100$false.neg.rate<-round(result.data1.exp5.100$fn/(result.data1.exp5.100$fn+result.data1.exp5.100$tp)*100,1)
result.data1.exp5.100$false.disc.rate<-round(result.data1.exp5.100$fp/(result.data1.exp5.100$fp+result.data1.exp5.100$tp)*100,1)
result.data1.exp5.100$false.omit.rate<-round(result.data1.exp5.100$fn/(result.data1.exp5.100$fn+result.data1.exp5.100$tn)*100,1)

result.data1.exp5.100$expo<-5
result.data1.exp5.100$ICC<-1


# Low ICC
result.data1.exp5.100.ICC.low<-NULL
for (j in 4:(dim(RES1.total.exp5.100.ICC.low)[2]-1)){
  RES1.total.exp5.100.ICC.low[,j]<-as.factor(as.character(RES1.total.exp5.100.ICC.low[,j]))  
  RES1.total.exp5.100.ICC.low$true.pred<-as.factor(as.character(RES1.total.exp5.100.ICC.low$true.pred))
  tn<-squareTable(RES1.total.exp5.100.ICC.low[,j],RES1.total.exp5.100.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp5.100.ICC.low[,j],RES1.total.exp5.100.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp5.100.ICC.low[,j],RES1.total.exp5.100.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp5.100.ICC.low[,j],RES1.total.exp5.100.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp5.100.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp5.100.ICC.low<-rbind(result.data1.exp5.100.ICC.low,tot)
}


result.data1.exp5.100.ICC.low<-as.data.frame(result.data1.exp5.100.ICC.low)
colnames(result.data1.exp5.100.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp5.100.ICC.low$nselected<-as.numeric(as.character(result.data1.exp5.100.ICC.low$nselected))
result.data1.exp5.100.ICC.low$tn<-as.numeric(as.character(result.data1.exp5.100.ICC.low$tn))
result.data1.exp5.100.ICC.low$fn<-as.numeric(as.character(result.data1.exp5.100.ICC.low$fn))
result.data1.exp5.100.ICC.low$fp<-as.numeric(as.character(result.data1.exp5.100.ICC.low$fp))
result.data1.exp5.100.ICC.low$tp<-as.numeric(as.character(result.data1.exp5.100.ICC.low$tp))

result.data1.exp5.100.ICC.low$sensitivity<-round(result.data1.exp5.100.ICC.low$tp/(result.data1.exp5.100.ICC.low$tp+result.data1.exp5.100.ICC.low$fn)*100,1)
result.data1.exp5.100.ICC.low$specificity<-round(result.data1.exp5.100.ICC.low$tn/(result.data1.exp5.100.ICC.low$tn+result.data1.exp5.100.ICC.low$fp)*100,1)
result.data1.exp5.100.ICC.low$false.pos.rate<-round(result.data1.exp5.100.ICC.low$fp/(result.data1.exp5.100.ICC.low$fp+result.data1.exp5.100.ICC.low$tn)*100,1)
result.data1.exp5.100.ICC.low$false.neg.rate<-round(result.data1.exp5.100.ICC.low$fn/(result.data1.exp5.100.ICC.low$fn+result.data1.exp5.100.ICC.low$tp)*100,1)
result.data1.exp5.100.ICC.low$false.disc.rate<-round(result.data1.exp5.100.ICC.low$fp/(result.data1.exp5.100.ICC.low$fp+result.data1.exp5.100.ICC.low$tp)*100,1)
result.data1.exp5.100.ICC.low$false.omit.rate<-round(result.data1.exp5.100.ICC.low$fn/(result.data1.exp5.100.ICC.low$fn+result.data1.exp5.100.ICC.low$tn)*100,1)

result.data1.exp5.100.ICC.low$expo<-5
result.data1.exp5.100.ICC.low$ICC<-0.1



# Medium ICC
result.data1.exp5.100.ICC.med<-NULL
for (j in 4:(dim(RES1.total.exp5.100.ICC.med)[2]-1)){
  RES1.total.exp5.100.ICC.med[,j]<-as.factor(as.character(RES1.total.exp5.100.ICC.med[,j]))  
  RES1.total.exp5.100.ICC.med$true.pred<-as.factor(as.character(RES1.total.exp5.100.ICC.med$true.pred))
  tn<-squareTable(RES1.total.exp5.100.ICC.med[,j],RES1.total.exp5.100.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp5.100.ICC.med[,j],RES1.total.exp5.100.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp5.100.ICC.med[,j],RES1.total.exp5.100.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp5.100.ICC.med[,j],RES1.total.exp5.100.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp5.100.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp5.100.ICC.med<-rbind(result.data1.exp5.100.ICC.med,tot)
}


result.data1.exp5.100.ICC.med<-as.data.frame(result.data1.exp5.100.ICC.med)
colnames(result.data1.exp5.100.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp5.100.ICC.med$nselected<-as.numeric(as.character(result.data1.exp5.100.ICC.med$nselected))
result.data1.exp5.100.ICC.med$tn<-as.numeric(as.character(result.data1.exp5.100.ICC.med$tn))
result.data1.exp5.100.ICC.med$fn<-as.numeric(as.character(result.data1.exp5.100.ICC.med$fn))
result.data1.exp5.100.ICC.med$fp<-as.numeric(as.character(result.data1.exp5.100.ICC.med$fp))
result.data1.exp5.100.ICC.med$tp<-as.numeric(as.character(result.data1.exp5.100.ICC.med$tp))

result.data1.exp5.100.ICC.med$sensitivity<-round(result.data1.exp5.100.ICC.med$tp/(result.data1.exp5.100.ICC.med$tp+result.data1.exp5.100.ICC.med$fn)*100,1)
result.data1.exp5.100.ICC.med$specificity<-round(result.data1.exp5.100.ICC.med$tn/(result.data1.exp5.100.ICC.med$tn+result.data1.exp5.100.ICC.med$fp)*100,1)
result.data1.exp5.100.ICC.med$false.pos.rate<-round(result.data1.exp5.100.ICC.med$fp/(result.data1.exp5.100.ICC.med$fp+result.data1.exp5.100.ICC.med$tn)*100,1)
result.data1.exp5.100.ICC.med$false.neg.rate<-round(result.data1.exp5.100.ICC.med$fn/(result.data1.exp5.100.ICC.med$fn+result.data1.exp5.100.ICC.med$tp)*100,1)
result.data1.exp5.100.ICC.med$false.disc.rate<-round(result.data1.exp5.100.ICC.med$fp/(result.data1.exp5.100.ICC.med$fp+result.data1.exp5.100.ICC.med$tp)*100,1)
result.data1.exp5.100.ICC.med$false.omit.rate<-round(result.data1.exp5.100.ICC.med$fn/(result.data1.exp5.100.ICC.med$fn+result.data1.exp5.100.ICC.med$tn)*100,1)

result.data1.exp5.100.ICC.med$expo<-5
result.data1.exp5.100.ICC.med$ICC<-0.5


# High ICC
result.data1.exp5.100.ICC.high<-NULL
for (j in 4:(dim(RES1.total.exp5.100.ICC.high)[2]-1)){
  RES1.total.exp5.100.ICC.high[,j]<-as.factor(as.character(RES1.total.exp5.100.ICC.high[,j]))  
  RES1.total.exp5.100.ICC.high$true.pred<-as.factor(as.character(RES1.total.exp5.100.ICC.high$true.pred))
  tn<-squareTable(RES1.total.exp5.100.ICC.high[,j],RES1.total.exp5.100.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp5.100.ICC.high[,j],RES1.total.exp5.100.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp5.100.ICC.high[,j],RES1.total.exp5.100.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp5.100.ICC.high[,j],RES1.total.exp5.100.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp5.100.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp5.100.ICC.high<-rbind(result.data1.exp5.100.ICC.high,tot)
}


result.data1.exp5.100.ICC.high<-as.data.frame(result.data1.exp5.100.ICC.high)
colnames(result.data1.exp5.100.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp5.100.ICC.high$nselected<-as.numeric(as.character(result.data1.exp5.100.ICC.high$nselected))
result.data1.exp5.100.ICC.high$tn<-as.numeric(as.character(result.data1.exp5.100.ICC.high$tn))
result.data1.exp5.100.ICC.high$fn<-as.numeric(as.character(result.data1.exp5.100.ICC.high$fn))
result.data1.exp5.100.ICC.high$fp<-as.numeric(as.character(result.data1.exp5.100.ICC.high$fp))
result.data1.exp5.100.ICC.high$tp<-as.numeric(as.character(result.data1.exp5.100.ICC.high$tp))

result.data1.exp5.100.ICC.high$sensitivity<-round(result.data1.exp5.100.ICC.high$tp/(result.data1.exp5.100.ICC.high$tp+result.data1.exp5.100.ICC.high$fn)*100,1)
result.data1.exp5.100.ICC.high$specificity<-round(result.data1.exp5.100.ICC.high$tn/(result.data1.exp5.100.ICC.high$tn+result.data1.exp5.100.ICC.high$fp)*100,1)
result.data1.exp5.100.ICC.high$false.pos.rate<-round(result.data1.exp5.100.ICC.high$fp/(result.data1.exp5.100.ICC.high$fp+result.data1.exp5.100.ICC.high$tn)*100,1)
result.data1.exp5.100.ICC.high$false.neg.rate<-round(result.data1.exp5.100.ICC.high$fn/(result.data1.exp5.100.ICC.high$fn+result.data1.exp5.100.ICC.high$tp)*100,1)
result.data1.exp5.100.ICC.high$false.disc.rate<-round(result.data1.exp5.100.ICC.high$fp/(result.data1.exp5.100.ICC.high$fp+result.data1.exp5.100.ICC.high$tp)*100,1)
result.data1.exp5.100.ICC.high$false.omit.rate<-round(result.data1.exp5.100.ICC.high$fn/(result.data1.exp5.100.ICC.high$fn+result.data1.exp5.100.ICC.high$tn)*100,1)

result.data1.exp5.100.ICC.high$expo<-5
result.data1.exp5.100.ICC.high$ICC<-0.9


############### exp10

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY1andX/exp10")
nexp<-10


load(file="RES1.all.exp10.100.2023.RData")
colnames(RES1.all.exp10.100)
colnames(RES1.all.exp10.100)<-c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES1.all.exp10.100<-RES1.all.exp10.100[,c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                                "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]

load(file="twostep/RES1avRed.all.exp10.100.RData")
colnames(RES1avRed.all.exp10.100)
colnames(RES1avRed.all.exp10.100)<-c("numsim","var","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                    "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES1.DLNM.all.exp10.100.RData")
colnames(RES1.DLNM.all.exp10.100)
colnames(RES1.DLNM.all.exp10.100)<-c("numsim","var","true.pred","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                    "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                    "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",     
                                    "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")


load(file="dlnm/RES1.DLNM.AVG.all.exp10.100.RData")
colnames(RES1.DLNM.AVG.all.exp10.100)
colnames(RES1.DLNM.AVG.all.exp10.100)<-c("numsim","var","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                        "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES1.total.exp10.100<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                            list(RES1.all.exp10.100,RES1avRed.all.exp10.100,RES1.DLNM.all.exp10.100,RES1.DLNM.AVG.all.exp10.100))


### Merge avec ICC

ICC.data1.exp10.red<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data1.exp10.red)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data1.exp10.red<-rbind(ICC.data1.exp10.red,ICC.i[[i]])
}
ICC.data1.exp10.red<-ICC.data1.exp10.red[-1,]
ICC.data1.exp10.red$numsim<-as.integer(ICC.data1.exp10.red$numsim)
colnames(ICC.data1.exp10.red)<-c("var","numsim","ICC")

ICC.data1.exp10.red$expo_name<-ifelse(nchar(as.character(ICC.data1.exp10.red$var))==4,substring(ICC.data1.exp10.red$var,first=1,last=2),
                                     ifelse(nchar(as.character(ICC.data1.exp10.red$var))==5,substring(ICC.data1.exp10.red$var,first=1,last=3),
                                            substring(ICC.data1.exp10.red$var,first=1,last=4)))

ICC.data1.exp10.red<-ICC.data1.exp10.red[,c("expo_name","numsim","ICC")]
ICC.data1.exp10.red<-unique(ICC.data1.exp10.red)
colnames(ICC.data1.exp10.red)<-c("var","numsim","ICC")

RES1.total.exp10.100.ICC<-merge(RES1.total.exp10.100,ICC.data1.exp10.red,by=c("var","numsim"))

table(RES1.total.exp10.100.ICC$true.pred,RES1.total.exp10.100.ICC$ICC)


squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


RES1.total.exp10.100.ICC.low<-RES1.total.exp10.100.ICC[RES1.total.exp10.100.ICC$ICC==0.1,]
RES1.total.exp10.100.ICC.med<-RES1.total.exp10.100.ICC[RES1.total.exp10.100.ICC$ICC==0.5,]
RES1.total.exp10.100.ICC.high<-RES1.total.exp10.100.ICC[RES1.total.exp10.100.ICC$ICC==0.9,]


### All
RES1.i.total.exp10.100.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES1.i.total.exp10.100.ICC[[i]]<-RES1.total.exp10.100.ICC[RES1.total.exp10.100.ICC$numsim == i,]
}

result.data1.i.exp10.100<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES1.i.total.exp10.100.ICC[[i]])[2]-1)){
    RES1.i.total.exp10.100.ICC[[i]][,j]<-as.factor(as.character(RES1.i.total.exp10.100.ICC[[i]][,j]))  
    RES1.i.total.exp10.100.ICC[[i]]$true.pred<-as.factor(as.character(RES1.i.total.exp10.100.ICC[[i]]$true.pred))
    tn<-squareTable(RES1.i.total.exp10.100.ICC[[i]][,j],RES1.i.total.exp10.100.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES1.i.total.exp10.100.ICC[[i]][,j],RES1.i.total.exp10.100.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES1.i.total.exp10.100.ICC[[i]][,j],RES1.i.total.exp10.100.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES1.i.total.exp10.100.ICC[[i]][,j],RES1.i.total.exp10.100.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES1.i.total.exp10.100.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data1.i.exp10.100[[i]]<-rbind(result.data1.i.exp10.100[[i]],tot)
  }
}


result.data1.exp10.100<-matrix(ncol=dim(result.data1.i.exp10.100[[1]])[2])
for (i in 1:nsim){
  result.data1.exp10.100<-rbind(result.data1.exp10.100,result.data1.i.exp10.100[[i]])
}

result.data1.exp10.100<-as.data.frame(result.data1.exp10.100)
colnames(result.data1.exp10.100)<-c("Method","nselected","tn","fn","fp","tp","data.i")


result.data1.exp10.100$nselected<-as.numeric(as.character(result.data1.exp10.100$nselected))
result.data1.exp10.100$tn<-as.numeric(as.character(result.data1.exp10.100$tn))
result.data1.exp10.100$fn<-as.numeric(as.character(result.data1.exp10.100$fn))
result.data1.exp10.100$fp<-as.numeric(as.character(result.data1.exp10.100$fp))
result.data1.exp10.100$tp<-as.numeric(as.character(result.data1.exp10.100$tp))

result.data1.exp10.100$sensitivity<-round(result.data1.exp10.100$tp/(result.data1.exp10.100$tp+result.data1.exp10.100$fn)*100,1)
result.data1.exp10.100$specificity<-round(result.data1.exp10.100$tn/(result.data1.exp10.100$tn+result.data1.exp10.100$fp)*100,1)
result.data1.exp10.100$false.pos.rate<-round(result.data1.exp10.100$fp/(result.data1.exp10.100$fp+result.data1.exp10.100$tn)*100,1)
result.data1.exp10.100$false.neg.rate<-round(result.data1.exp10.100$fn/(result.data1.exp10.100$fn+result.data1.exp10.100$tp)*100,1)
result.data1.exp10.100$false.disc.rate<-round(result.data1.exp10.100$fp/(result.data1.exp10.100$fp+result.data1.exp10.100$tp)*100,1)
result.data1.exp10.100$false.omit.rate<-round(result.data1.exp10.100$fn/(result.data1.exp10.100$fn+result.data1.exp10.100$tn)*100,1)

result.data1.exp10.100$expo<-10
result.data1.exp10.100$ICC<-1


# Low ICC
result.data1.exp10.100.ICC.low<-NULL
for (j in 4:(dim(RES1.total.exp10.100.ICC.low)[2]-1)){
  RES1.total.exp10.100.ICC.low[,j]<-as.factor(as.character(RES1.total.exp10.100.ICC.low[,j]))  
  RES1.total.exp10.100.ICC.low$true.pred<-as.factor(as.character(RES1.total.exp10.100.ICC.low$true.pred))
  tn<-squareTable(RES1.total.exp10.100.ICC.low[,j],RES1.total.exp10.100.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp10.100.ICC.low[,j],RES1.total.exp10.100.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp10.100.ICC.low[,j],RES1.total.exp10.100.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp10.100.ICC.low[,j],RES1.total.exp10.100.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp10.100.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp10.100.ICC.low<-rbind(result.data1.exp10.100.ICC.low,tot)
}


result.data1.exp10.100.ICC.low<-as.data.frame(result.data1.exp10.100.ICC.low)
colnames(result.data1.exp10.100.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp10.100.ICC.low$nselected<-as.numeric(as.character(result.data1.exp10.100.ICC.low$nselected))
result.data1.exp10.100.ICC.low$tn<-as.numeric(as.character(result.data1.exp10.100.ICC.low$tn))
result.data1.exp10.100.ICC.low$fn<-as.numeric(as.character(result.data1.exp10.100.ICC.low$fn))
result.data1.exp10.100.ICC.low$fp<-as.numeric(as.character(result.data1.exp10.100.ICC.low$fp))
result.data1.exp10.100.ICC.low$tp<-as.numeric(as.character(result.data1.exp10.100.ICC.low$tp))

result.data1.exp10.100.ICC.low$sensitivity<-round(result.data1.exp10.100.ICC.low$tp/(result.data1.exp10.100.ICC.low$tp+result.data1.exp10.100.ICC.low$fn)*100,1)
result.data1.exp10.100.ICC.low$specificity<-round(result.data1.exp10.100.ICC.low$tn/(result.data1.exp10.100.ICC.low$tn+result.data1.exp10.100.ICC.low$fp)*100,1)
result.data1.exp10.100.ICC.low$false.pos.rate<-round(result.data1.exp10.100.ICC.low$fp/(result.data1.exp10.100.ICC.low$fp+result.data1.exp10.100.ICC.low$tn)*100,1)
result.data1.exp10.100.ICC.low$false.neg.rate<-round(result.data1.exp10.100.ICC.low$fn/(result.data1.exp10.100.ICC.low$fn+result.data1.exp10.100.ICC.low$tp)*100,1)
result.data1.exp10.100.ICC.low$false.disc.rate<-round(result.data1.exp10.100.ICC.low$fp/(result.data1.exp10.100.ICC.low$fp+result.data1.exp10.100.ICC.low$tp)*100,1)
result.data1.exp10.100.ICC.low$false.omit.rate<-round(result.data1.exp10.100.ICC.low$fn/(result.data1.exp10.100.ICC.low$fn+result.data1.exp10.100.ICC.low$tn)*100,1)

result.data1.exp10.100.ICC.low$expo<-10
result.data1.exp10.100.ICC.low$ICC<-0.1



# Medium ICC
result.data1.exp10.100.ICC.med<-NULL
for (j in 4:(dim(RES1.total.exp10.100.ICC.med)[2]-1)){
  RES1.total.exp10.100.ICC.med[,j]<-as.factor(as.character(RES1.total.exp10.100.ICC.med[,j]))  
  RES1.total.exp10.100.ICC.med$true.pred<-as.factor(as.character(RES1.total.exp10.100.ICC.med$true.pred))
  tn<-squareTable(RES1.total.exp10.100.ICC.med[,j],RES1.total.exp10.100.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp10.100.ICC.med[,j],RES1.total.exp10.100.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp10.100.ICC.med[,j],RES1.total.exp10.100.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp10.100.ICC.med[,j],RES1.total.exp10.100.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp10.100.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp10.100.ICC.med<-rbind(result.data1.exp10.100.ICC.med,tot)
}


result.data1.exp10.100.ICC.med<-as.data.frame(result.data1.exp10.100.ICC.med)
colnames(result.data1.exp10.100.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp10.100.ICC.med$nselected<-as.numeric(as.character(result.data1.exp10.100.ICC.med$nselected))
result.data1.exp10.100.ICC.med$tn<-as.numeric(as.character(result.data1.exp10.100.ICC.med$tn))
result.data1.exp10.100.ICC.med$fn<-as.numeric(as.character(result.data1.exp10.100.ICC.med$fn))
result.data1.exp10.100.ICC.med$fp<-as.numeric(as.character(result.data1.exp10.100.ICC.med$fp))
result.data1.exp10.100.ICC.med$tp<-as.numeric(as.character(result.data1.exp10.100.ICC.med$tp))

result.data1.exp10.100.ICC.med$sensitivity<-round(result.data1.exp10.100.ICC.med$tp/(result.data1.exp10.100.ICC.med$tp+result.data1.exp10.100.ICC.med$fn)*100,1)
result.data1.exp10.100.ICC.med$specificity<-round(result.data1.exp10.100.ICC.med$tn/(result.data1.exp10.100.ICC.med$tn+result.data1.exp10.100.ICC.med$fp)*100,1)
result.data1.exp10.100.ICC.med$false.pos.rate<-round(result.data1.exp10.100.ICC.med$fp/(result.data1.exp10.100.ICC.med$fp+result.data1.exp10.100.ICC.med$tn)*100,1)
result.data1.exp10.100.ICC.med$false.neg.rate<-round(result.data1.exp10.100.ICC.med$fn/(result.data1.exp10.100.ICC.med$fn+result.data1.exp10.100.ICC.med$tp)*100,1)
result.data1.exp10.100.ICC.med$false.disc.rate<-round(result.data1.exp10.100.ICC.med$fp/(result.data1.exp10.100.ICC.med$fp+result.data1.exp10.100.ICC.med$tp)*100,1)
result.data1.exp10.100.ICC.med$false.omit.rate<-round(result.data1.exp10.100.ICC.med$fn/(result.data1.exp10.100.ICC.med$fn+result.data1.exp10.100.ICC.med$tn)*100,1)

result.data1.exp10.100.ICC.med$expo<-10
result.data1.exp10.100.ICC.med$ICC<-0.5


# High ICC
result.data1.exp10.100.ICC.high<-NULL
for (j in 4:(dim(RES1.total.exp10.100.ICC.high)[2]-1)){
  RES1.total.exp10.100.ICC.high[,j]<-as.factor(as.character(RES1.total.exp10.100.ICC.high[,j]))  
  RES1.total.exp10.100.ICC.high$true.pred<-as.factor(as.character(RES1.total.exp10.100.ICC.high$true.pred))
  tn<-squareTable(RES1.total.exp10.100.ICC.high[,j],RES1.total.exp10.100.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES1.total.exp10.100.ICC.high[,j],RES1.total.exp10.100.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES1.total.exp10.100.ICC.high[,j],RES1.total.exp10.100.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES1.total.exp10.100.ICC.high[,j],RES1.total.exp10.100.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES1.total.exp10.100.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data1.exp10.100.ICC.high<-rbind(result.data1.exp10.100.ICC.high,tot)
}


result.data1.exp10.100.ICC.high<-as.data.frame(result.data1.exp10.100.ICC.high)
colnames(result.data1.exp10.100.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data1.exp10.100.ICC.high$nselected<-as.numeric(as.character(result.data1.exp10.100.ICC.high$nselected))
result.data1.exp10.100.ICC.high$tn<-as.numeric(as.character(result.data1.exp10.100.ICC.high$tn))
result.data1.exp10.100.ICC.high$fn<-as.numeric(as.character(result.data1.exp10.100.ICC.high$fn))
result.data1.exp10.100.ICC.high$fp<-as.numeric(as.character(result.data1.exp10.100.ICC.high$fp))
result.data1.exp10.100.ICC.high$tp<-as.numeric(as.character(result.data1.exp10.100.ICC.high$tp))

result.data1.exp10.100.ICC.high$sensitivity<-round(result.data1.exp10.100.ICC.high$tp/(result.data1.exp10.100.ICC.high$tp+result.data1.exp10.100.ICC.high$fn)*100,1)
result.data1.exp10.100.ICC.high$specificity<-round(result.data1.exp10.100.ICC.high$tn/(result.data1.exp10.100.ICC.high$tn+result.data1.exp10.100.ICC.high$fp)*100,1)
result.data1.exp10.100.ICC.high$false.pos.rate<-round(result.data1.exp10.100.ICC.high$fp/(result.data1.exp10.100.ICC.high$fp+result.data1.exp10.100.ICC.high$tn)*100,1)
result.data1.exp10.100.ICC.high$false.neg.rate<-round(result.data1.exp10.100.ICC.high$fn/(result.data1.exp10.100.ICC.high$fn+result.data1.exp10.100.ICC.high$tp)*100,1)
result.data1.exp10.100.ICC.high$false.disc.rate<-round(result.data1.exp10.100.ICC.high$fp/(result.data1.exp10.100.ICC.high$fp+result.data1.exp10.100.ICC.high$tp)*100,1)
result.data1.exp10.100.ICC.high$false.omit.rate<-round(result.data1.exp10.100.ICC.high$fn/(result.data1.exp10.100.ICC.high$fn+result.data1.exp10.100.ICC.high$tn)*100,1)

result.data1.exp10.100.ICC.high$expo<-10
result.data1.exp10.100.ICC.high$ICC<-0.9


result.data1.100.ICC<-rbind(result.data1.exp3.100,result.data1.exp3.100.ICC.low,result.data1.exp3.100.ICC.med,result.data1.exp3.100.ICC.high,
                            result.data1.exp5.100,result.data1.exp5.100.ICC.low,result.data1.exp5.100.ICC.med,result.data1.exp5.100.ICC.high,
                            result.data1.exp10.100,result.data1.exp10.100.ICC.low,result.data1.exp10.100.ICC.med,result.data1.exp10.100.ICC.high)

result.data1.100.ICC$false.disc.rate[result.data1.100.ICC$nselected==0]<-0

save(result.data1.100.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data1.100_ICC.2023.Rdata")
write.csv2(result.data1.100.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary.100_data1_ICC.2023.csv")




############################################# data2

############### EXP3

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY2andX/exp3")
nexp<-3


load(file="RES2.all.exp3.100.2023.RData")
colnames(RES2.all.exp3.100)
colnames(RES2.all.exp3.100)<-c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES2.all.exp3.100<-RES2.all.exp3.100[,c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES2avRed.all.exp3.100.RData")
colnames(RES2avRed.all.exp3.100)
colnames(RES2avRed.all.exp3.100)<-c("numsim","var","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                    "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES2.DLNM.all.exp3.100.RData")
colnames(RES2.DLNM.all.exp3.100)
colnames(RES2.DLNM.all.exp3.100)<-c("numsim","var","true.pred","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                    "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                    "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",     
                                    "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")


load(file="dlnm/RES2.DLNM.AVG.all.exp3.100.RData")
colnames(RES2.DLNM.AVG.all.exp3.100)
colnames(RES2.DLNM.AVG.all.exp3.100)<-c("numsim","var","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                        "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES2.total.exp3.100<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                            list(RES2.all.exp3.100,RES2avRed.all.exp3.100,RES2.DLNM.all.exp3.100,RES2.DLNM.AVG.all.exp3.100))


### Merge avec ICC

ICC.data2.exp3.red<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data2.exp3.red)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data2.exp3.red<-rbind(ICC.data2.exp3.red,ICC.i[[i]])
}
ICC.data2.exp3.red<-ICC.data2.exp3.red[-1,]
ICC.data2.exp3.red$numsim<-as.integer(ICC.data2.exp3.red$numsim)
colnames(ICC.data2.exp3.red)<-c("var","numsim","ICC")

ICC.data2.exp3.red$expo_name<-ifelse(nchar(as.character(ICC.data2.exp3.red$var))==4,substring(ICC.data2.exp3.red$var,first=1,last=2),
                                     ifelse(nchar(as.character(ICC.data2.exp3.red$var))==5,substring(ICC.data2.exp3.red$var,first=1,last=3),
                                            substring(ICC.data2.exp3.red$var,first=1,last=4)))

ICC.data2.exp3.red<-ICC.data2.exp3.red[,c("expo_name","numsim","ICC")]
ICC.data2.exp3.red<-unique(ICC.data2.exp3.red)
colnames(ICC.data2.exp3.red)<-c("var","numsim","ICC")

RES2.total.exp3.100.ICC<-merge(RES2.total.exp3.100,ICC.data2.exp3.red,by=c("var","numsim"))

table(RES2.total.exp3.100.ICC$true.pred,RES2.total.exp3.100.ICC$ICC)


squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


RES2.total.exp3.100.ICC.low<-RES2.total.exp3.100.ICC[RES2.total.exp3.100.ICC$ICC==0.1,]
RES2.total.exp3.100.ICC.med<-RES2.total.exp3.100.ICC[RES2.total.exp3.100.ICC$ICC==0.5,]
RES2.total.exp3.100.ICC.high<-RES2.total.exp3.100.ICC[RES2.total.exp3.100.ICC$ICC==0.9,]


### All
RES2.i.total.exp3.100.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES2.i.total.exp3.100.ICC[[i]]<-RES2.total.exp3.100.ICC[RES2.total.exp3.100.ICC$numsim == i,]
}

result.data2.i.exp3.100<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES2.i.total.exp3.100.ICC[[i]])[2]-1)){
    RES2.i.total.exp3.100.ICC[[i]][,j]<-as.factor(as.character(RES2.i.total.exp3.100.ICC[[i]][,j]))  
    RES2.i.total.exp3.100.ICC[[i]]$true.pred<-as.factor(as.character(RES2.i.total.exp3.100.ICC[[i]]$true.pred))
    tn<-squareTable(RES2.i.total.exp3.100.ICC[[i]][,j],RES2.i.total.exp3.100.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES2.i.total.exp3.100.ICC[[i]][,j],RES2.i.total.exp3.100.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES2.i.total.exp3.100.ICC[[i]][,j],RES2.i.total.exp3.100.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES2.i.total.exp3.100.ICC[[i]][,j],RES2.i.total.exp3.100.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES2.i.total.exp3.100.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data2.i.exp3.100[[i]]<-rbind(result.data2.i.exp3.100[[i]],tot)
  }
}


result.data2.exp3.100<-matrix(ncol=dim(result.data2.i.exp3.100[[1]])[2])
for (i in 1:nsim){
  result.data2.exp3.100<-rbind(result.data2.exp3.100,result.data2.i.exp3.100[[i]])
}

result.data2.exp3.100<-as.data.frame(result.data2.exp3.100)
colnames(result.data2.exp3.100)<-c("Method","nselected","tn","fn","fp","tp","data.i")


result.data2.exp3.100$nselected<-as.numeric(as.character(result.data2.exp3.100$nselected))
result.data2.exp3.100$tn<-as.numeric(as.character(result.data2.exp3.100$tn))
result.data2.exp3.100$fn<-as.numeric(as.character(result.data2.exp3.100$fn))
result.data2.exp3.100$fp<-as.numeric(as.character(result.data2.exp3.100$fp))
result.data2.exp3.100$tp<-as.numeric(as.character(result.data2.exp3.100$tp))

result.data2.exp3.100$sensitivity<-round(result.data2.exp3.100$tp/(result.data2.exp3.100$tp+result.data2.exp3.100$fn)*100,1)
result.data2.exp3.100$specificity<-round(result.data2.exp3.100$tn/(result.data2.exp3.100$tn+result.data2.exp3.100$fp)*100,1)
result.data2.exp3.100$false.pos.rate<-round(result.data2.exp3.100$fp/(result.data2.exp3.100$fp+result.data2.exp3.100$tn)*100,1)
result.data2.exp3.100$false.neg.rate<-round(result.data2.exp3.100$fn/(result.data2.exp3.100$fn+result.data2.exp3.100$tp)*100,1)
result.data2.exp3.100$false.disc.rate<-round(result.data2.exp3.100$fp/(result.data2.exp3.100$fp+result.data2.exp3.100$tp)*100,1)
result.data2.exp3.100$false.omit.rate<-round(result.data2.exp3.100$fn/(result.data2.exp3.100$fn+result.data2.exp3.100$tn)*100,1)

result.data2.exp3.100$expo<-3
result.data2.exp3.100$ICC<-1


# Low ICC
result.data2.exp3.100.ICC.low<-NULL
for (j in 4:(dim(RES2.total.exp3.100.ICC.low)[2]-1)){
  RES2.total.exp3.100.ICC.low[,j]<-as.factor(as.character(RES2.total.exp3.100.ICC.low[,j]))  
  RES2.total.exp3.100.ICC.low$true.pred<-as.factor(as.character(RES2.total.exp3.100.ICC.low$true.pred))
  tn<-squareTable(RES2.total.exp3.100.ICC.low[,j],RES2.total.exp3.100.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp3.100.ICC.low[,j],RES2.total.exp3.100.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp3.100.ICC.low[,j],RES2.total.exp3.100.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp3.100.ICC.low[,j],RES2.total.exp3.100.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp3.100.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp3.100.ICC.low<-rbind(result.data2.exp3.100.ICC.low,tot)
}


result.data2.exp3.100.ICC.low<-as.data.frame(result.data2.exp3.100.ICC.low)
colnames(result.data2.exp3.100.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp3.100.ICC.low$nselected<-as.numeric(as.character(result.data2.exp3.100.ICC.low$nselected))
result.data2.exp3.100.ICC.low$tn<-as.numeric(as.character(result.data2.exp3.100.ICC.low$tn))
result.data2.exp3.100.ICC.low$fn<-as.numeric(as.character(result.data2.exp3.100.ICC.low$fn))
result.data2.exp3.100.ICC.low$fp<-as.numeric(as.character(result.data2.exp3.100.ICC.low$fp))
result.data2.exp3.100.ICC.low$tp<-as.numeric(as.character(result.data2.exp3.100.ICC.low$tp))

result.data2.exp3.100.ICC.low$sensitivity<-round(result.data2.exp3.100.ICC.low$tp/(result.data2.exp3.100.ICC.low$tp+result.data2.exp3.100.ICC.low$fn)*100,1)
result.data2.exp3.100.ICC.low$specificity<-round(result.data2.exp3.100.ICC.low$tn/(result.data2.exp3.100.ICC.low$tn+result.data2.exp3.100.ICC.low$fp)*100,1)
result.data2.exp3.100.ICC.low$false.pos.rate<-round(result.data2.exp3.100.ICC.low$fp/(result.data2.exp3.100.ICC.low$fp+result.data2.exp3.100.ICC.low$tn)*100,1)
result.data2.exp3.100.ICC.low$false.neg.rate<-round(result.data2.exp3.100.ICC.low$fn/(result.data2.exp3.100.ICC.low$fn+result.data2.exp3.100.ICC.low$tp)*100,1)
result.data2.exp3.100.ICC.low$false.disc.rate<-round(result.data2.exp3.100.ICC.low$fp/(result.data2.exp3.100.ICC.low$fp+result.data2.exp3.100.ICC.low$tp)*100,1)
result.data2.exp3.100.ICC.low$false.omit.rate<-round(result.data2.exp3.100.ICC.low$fn/(result.data2.exp3.100.ICC.low$fn+result.data2.exp3.100.ICC.low$tn)*100,1)

result.data2.exp3.100.ICC.low$expo<-3
result.data2.exp3.100.ICC.low$ICC<-0.1



# Medium ICC
result.data2.exp3.100.ICC.med<-NULL
for (j in 4:(dim(RES2.total.exp3.100.ICC.med)[2]-1)){
  RES2.total.exp3.100.ICC.med[,j]<-as.factor(as.character(RES2.total.exp3.100.ICC.med[,j]))  
  RES2.total.exp3.100.ICC.med$true.pred<-as.factor(as.character(RES2.total.exp3.100.ICC.med$true.pred))
  tn<-squareTable(RES2.total.exp3.100.ICC.med[,j],RES2.total.exp3.100.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp3.100.ICC.med[,j],RES2.total.exp3.100.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp3.100.ICC.med[,j],RES2.total.exp3.100.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp3.100.ICC.med[,j],RES2.total.exp3.100.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp3.100.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp3.100.ICC.med<-rbind(result.data2.exp3.100.ICC.med,tot)
}


result.data2.exp3.100.ICC.med<-as.data.frame(result.data2.exp3.100.ICC.med)
colnames(result.data2.exp3.100.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp3.100.ICC.med$nselected<-as.numeric(as.character(result.data2.exp3.100.ICC.med$nselected))
result.data2.exp3.100.ICC.med$tn<-as.numeric(as.character(result.data2.exp3.100.ICC.med$tn))
result.data2.exp3.100.ICC.med$fn<-as.numeric(as.character(result.data2.exp3.100.ICC.med$fn))
result.data2.exp3.100.ICC.med$fp<-as.numeric(as.character(result.data2.exp3.100.ICC.med$fp))
result.data2.exp3.100.ICC.med$tp<-as.numeric(as.character(result.data2.exp3.100.ICC.med$tp))

result.data2.exp3.100.ICC.med$sensitivity<-round(result.data2.exp3.100.ICC.med$tp/(result.data2.exp3.100.ICC.med$tp+result.data2.exp3.100.ICC.med$fn)*100,1)
result.data2.exp3.100.ICC.med$specificity<-round(result.data2.exp3.100.ICC.med$tn/(result.data2.exp3.100.ICC.med$tn+result.data2.exp3.100.ICC.med$fp)*100,1)
result.data2.exp3.100.ICC.med$false.pos.rate<-round(result.data2.exp3.100.ICC.med$fp/(result.data2.exp3.100.ICC.med$fp+result.data2.exp3.100.ICC.med$tn)*100,1)
result.data2.exp3.100.ICC.med$false.neg.rate<-round(result.data2.exp3.100.ICC.med$fn/(result.data2.exp3.100.ICC.med$fn+result.data2.exp3.100.ICC.med$tp)*100,1)
result.data2.exp3.100.ICC.med$false.disc.rate<-round(result.data2.exp3.100.ICC.med$fp/(result.data2.exp3.100.ICC.med$fp+result.data2.exp3.100.ICC.med$tp)*100,1)
result.data2.exp3.100.ICC.med$false.omit.rate<-round(result.data2.exp3.100.ICC.med$fn/(result.data2.exp3.100.ICC.med$fn+result.data2.exp3.100.ICC.med$tn)*100,1)

result.data2.exp3.100.ICC.med$expo<-3
result.data2.exp3.100.ICC.med$ICC<-0.5


# High ICC
result.data2.exp3.100.ICC.high<-NULL
for (j in 4:(dim(RES2.total.exp3.100.ICC.high)[2]-1)){
  RES2.total.exp3.100.ICC.high[,j]<-as.factor(as.character(RES2.total.exp3.100.ICC.high[,j]))  
  RES2.total.exp3.100.ICC.high$true.pred<-as.factor(as.character(RES2.total.exp3.100.ICC.high$true.pred))
  tn<-squareTable(RES2.total.exp3.100.ICC.high[,j],RES2.total.exp3.100.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp3.100.ICC.high[,j],RES2.total.exp3.100.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp3.100.ICC.high[,j],RES2.total.exp3.100.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp3.100.ICC.high[,j],RES2.total.exp3.100.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp3.100.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp3.100.ICC.high<-rbind(result.data2.exp3.100.ICC.high,tot)
}


result.data2.exp3.100.ICC.high<-as.data.frame(result.data2.exp3.100.ICC.high)
colnames(result.data2.exp3.100.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp3.100.ICC.high$nselected<-as.numeric(as.character(result.data2.exp3.100.ICC.high$nselected))
result.data2.exp3.100.ICC.high$tn<-as.numeric(as.character(result.data2.exp3.100.ICC.high$tn))
result.data2.exp3.100.ICC.high$fn<-as.numeric(as.character(result.data2.exp3.100.ICC.high$fn))
result.data2.exp3.100.ICC.high$fp<-as.numeric(as.character(result.data2.exp3.100.ICC.high$fp))
result.data2.exp3.100.ICC.high$tp<-as.numeric(as.character(result.data2.exp3.100.ICC.high$tp))

result.data2.exp3.100.ICC.high$sensitivity<-round(result.data2.exp3.100.ICC.high$tp/(result.data2.exp3.100.ICC.high$tp+result.data2.exp3.100.ICC.high$fn)*100,1)
result.data2.exp3.100.ICC.high$specificity<-round(result.data2.exp3.100.ICC.high$tn/(result.data2.exp3.100.ICC.high$tn+result.data2.exp3.100.ICC.high$fp)*100,1)
result.data2.exp3.100.ICC.high$false.pos.rate<-round(result.data2.exp3.100.ICC.high$fp/(result.data2.exp3.100.ICC.high$fp+result.data2.exp3.100.ICC.high$tn)*100,1)
result.data2.exp3.100.ICC.high$false.neg.rate<-round(result.data2.exp3.100.ICC.high$fn/(result.data2.exp3.100.ICC.high$fn+result.data2.exp3.100.ICC.high$tp)*100,1)
result.data2.exp3.100.ICC.high$false.disc.rate<-round(result.data2.exp3.100.ICC.high$fp/(result.data2.exp3.100.ICC.high$fp+result.data2.exp3.100.ICC.high$tp)*100,1)
result.data2.exp3.100.ICC.high$false.omit.rate<-round(result.data2.exp3.100.ICC.high$fn/(result.data2.exp3.100.ICC.high$fn+result.data2.exp3.100.ICC.high$tn)*100,1)

result.data2.exp3.100.ICC.high$expo<-3
result.data2.exp3.100.ICC.high$ICC<-0.9


############### exp5

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY2andX/exp5")
nexp<-5


load(file="RES2.all.exp5.100.2023.RData")
colnames(RES2.all.exp5.100)
colnames(RES2.all.exp5.100)<-c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES2.all.exp5.100<-RES2.all.exp5.100[,c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                               "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES2avRed.all.exp5.100.RData")
colnames(RES2avRed.all.exp5.100)
colnames(RES2avRed.all.exp5.100)<-c("numsim","var","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                    "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES2.DLNM.all.exp5.100.RData")
colnames(RES2.DLNM.all.exp5.100)
colnames(RES2.DLNM.all.exp5.100)<-c("numsim","var","true.pred","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                    "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                    "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",     
                                    "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")


load(file="dlnm/RES2.DLNM.AVG.all.exp5.100.RData")
colnames(RES2.DLNM.AVG.all.exp5.100)
colnames(RES2.DLNM.AVG.all.exp5.100)<-c("numsim","var","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                        "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES2.total.exp5.100<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                            list(RES2.all.exp5.100,RES2avRed.all.exp5.100,RES2.DLNM.all.exp5.100,RES2.DLNM.AVG.all.exp5.100))


### Merge avec ICC

ICC.data2.exp5.red<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data2.exp5.red)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data2.exp5.red<-rbind(ICC.data2.exp5.red,ICC.i[[i]])
}
ICC.data2.exp5.red<-ICC.data2.exp5.red[-1,]
ICC.data2.exp5.red$numsim<-as.integer(ICC.data2.exp5.red$numsim)
colnames(ICC.data2.exp5.red)<-c("var","numsim","ICC")

ICC.data2.exp5.red$expo_name<-ifelse(nchar(as.character(ICC.data2.exp5.red$var))==4,substring(ICC.data2.exp5.red$var,first=1,last=2),
                                     ifelse(nchar(as.character(ICC.data2.exp5.red$var))==5,substring(ICC.data2.exp5.red$var,first=1,last=3),
                                            substring(ICC.data2.exp5.red$var,first=1,last=4)))

ICC.data2.exp5.red<-ICC.data2.exp5.red[,c("expo_name","numsim","ICC")]
ICC.data2.exp5.red<-unique(ICC.data2.exp5.red)
colnames(ICC.data2.exp5.red)<-c("var","numsim","ICC")

RES2.total.exp5.100.ICC<-merge(RES2.total.exp5.100,ICC.data2.exp5.red,by=c("var","numsim"))

table(RES2.total.exp5.100.ICC$true.pred,RES2.total.exp5.100.ICC$ICC)


squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


RES2.total.exp5.100.ICC.low<-RES2.total.exp5.100.ICC[RES2.total.exp5.100.ICC$ICC==0.1,]
RES2.total.exp5.100.ICC.med<-RES2.total.exp5.100.ICC[RES2.total.exp5.100.ICC$ICC==0.5,]
RES2.total.exp5.100.ICC.high<-RES2.total.exp5.100.ICC[RES2.total.exp5.100.ICC$ICC==0.9,]


### All
RES2.i.total.exp5.100.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES2.i.total.exp5.100.ICC[[i]]<-RES2.total.exp5.100.ICC[RES2.total.exp5.100.ICC$numsim == i,]
}

result.data2.i.exp5.100<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES2.i.total.exp5.100.ICC[[i]])[2]-1)){
    RES2.i.total.exp5.100.ICC[[i]][,j]<-as.factor(as.character(RES2.i.total.exp5.100.ICC[[i]][,j]))  
    RES2.i.total.exp5.100.ICC[[i]]$true.pred<-as.factor(as.character(RES2.i.total.exp5.100.ICC[[i]]$true.pred))
    tn<-squareTable(RES2.i.total.exp5.100.ICC[[i]][,j],RES2.i.total.exp5.100.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES2.i.total.exp5.100.ICC[[i]][,j],RES2.i.total.exp5.100.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES2.i.total.exp5.100.ICC[[i]][,j],RES2.i.total.exp5.100.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES2.i.total.exp5.100.ICC[[i]][,j],RES2.i.total.exp5.100.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES2.i.total.exp5.100.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data2.i.exp5.100[[i]]<-rbind(result.data2.i.exp5.100[[i]],tot)
  }
}


result.data2.exp5.100<-matrix(ncol=dim(result.data2.i.exp5.100[[1]])[2])
for (i in 1:nsim){
  result.data2.exp5.100<-rbind(result.data2.exp5.100,result.data2.i.exp5.100[[i]])
}

result.data2.exp5.100<-as.data.frame(result.data2.exp5.100)
colnames(result.data2.exp5.100)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.100$nselected<-as.numeric(as.character(result.data2.exp5.100$nselected))
result.data2.exp5.100$tn<-as.numeric(as.character(result.data2.exp5.100$tn))
result.data2.exp5.100$fn<-as.numeric(as.character(result.data2.exp5.100$fn))
result.data2.exp5.100$fp<-as.numeric(as.character(result.data2.exp5.100$fp))
result.data2.exp5.100$tp<-as.numeric(as.character(result.data2.exp5.100$tp))

result.data2.exp5.100$sensitivity<-round(result.data2.exp5.100$tp/(result.data2.exp5.100$tp+result.data2.exp5.100$fn)*100,1)
result.data2.exp5.100$specificity<-round(result.data2.exp5.100$tn/(result.data2.exp5.100$tn+result.data2.exp5.100$fp)*100,1)
result.data2.exp5.100$false.pos.rate<-round(result.data2.exp5.100$fp/(result.data2.exp5.100$fp+result.data2.exp5.100$tn)*100,1)
result.data2.exp5.100$false.neg.rate<-round(result.data2.exp5.100$fn/(result.data2.exp5.100$fn+result.data2.exp5.100$tp)*100,1)
result.data2.exp5.100$false.disc.rate<-round(result.data2.exp5.100$fp/(result.data2.exp5.100$fp+result.data2.exp5.100$tp)*100,1)
result.data2.exp5.100$false.omit.rate<-round(result.data2.exp5.100$fn/(result.data2.exp5.100$fn+result.data2.exp5.100$tn)*100,1)

result.data2.exp5.100$expo<-5
result.data2.exp5.100$ICC<-1


# Low ICC
result.data2.exp5.100.ICC.low<-NULL
for (j in 4:(dim(RES2.total.exp5.100.ICC.low)[2]-1)){
  RES2.total.exp5.100.ICC.low[,j]<-as.factor(as.character(RES2.total.exp5.100.ICC.low[,j]))  
  RES2.total.exp5.100.ICC.low$true.pred<-as.factor(as.character(RES2.total.exp5.100.ICC.low$true.pred))
  tn<-squareTable(RES2.total.exp5.100.ICC.low[,j],RES2.total.exp5.100.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp5.100.ICC.low[,j],RES2.total.exp5.100.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp5.100.ICC.low[,j],RES2.total.exp5.100.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp5.100.ICC.low[,j],RES2.total.exp5.100.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp5.100.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp5.100.ICC.low<-rbind(result.data2.exp5.100.ICC.low,tot)
}


result.data2.exp5.100.ICC.low<-as.data.frame(result.data2.exp5.100.ICC.low)
colnames(result.data2.exp5.100.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.100.ICC.low$nselected<-as.numeric(as.character(result.data2.exp5.100.ICC.low$nselected))
result.data2.exp5.100.ICC.low$tn<-as.numeric(as.character(result.data2.exp5.100.ICC.low$tn))
result.data2.exp5.100.ICC.low$fn<-as.numeric(as.character(result.data2.exp5.100.ICC.low$fn))
result.data2.exp5.100.ICC.low$fp<-as.numeric(as.character(result.data2.exp5.100.ICC.low$fp))
result.data2.exp5.100.ICC.low$tp<-as.numeric(as.character(result.data2.exp5.100.ICC.low$tp))

result.data2.exp5.100.ICC.low$sensitivity<-round(result.data2.exp5.100.ICC.low$tp/(result.data2.exp5.100.ICC.low$tp+result.data2.exp5.100.ICC.low$fn)*100,1)
result.data2.exp5.100.ICC.low$specificity<-round(result.data2.exp5.100.ICC.low$tn/(result.data2.exp5.100.ICC.low$tn+result.data2.exp5.100.ICC.low$fp)*100,1)
result.data2.exp5.100.ICC.low$false.pos.rate<-round(result.data2.exp5.100.ICC.low$fp/(result.data2.exp5.100.ICC.low$fp+result.data2.exp5.100.ICC.low$tn)*100,1)
result.data2.exp5.100.ICC.low$false.neg.rate<-round(result.data2.exp5.100.ICC.low$fn/(result.data2.exp5.100.ICC.low$fn+result.data2.exp5.100.ICC.low$tp)*100,1)
result.data2.exp5.100.ICC.low$false.disc.rate<-round(result.data2.exp5.100.ICC.low$fp/(result.data2.exp5.100.ICC.low$fp+result.data2.exp5.100.ICC.low$tp)*100,1)
result.data2.exp5.100.ICC.low$false.omit.rate<-round(result.data2.exp5.100.ICC.low$fn/(result.data2.exp5.100.ICC.low$fn+result.data2.exp5.100.ICC.low$tn)*100,1)

result.data2.exp5.100.ICC.low$expo<-5
result.data2.exp5.100.ICC.low$ICC<-0.1



# Medium ICC
result.data2.exp5.100.ICC.med<-NULL
for (j in 4:(dim(RES2.total.exp5.100.ICC.med)[2]-1)){
  RES2.total.exp5.100.ICC.med[,j]<-as.factor(as.character(RES2.total.exp5.100.ICC.med[,j]))  
  RES2.total.exp5.100.ICC.med$true.pred<-as.factor(as.character(RES2.total.exp5.100.ICC.med$true.pred))
  tn<-squareTable(RES2.total.exp5.100.ICC.med[,j],RES2.total.exp5.100.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp5.100.ICC.med[,j],RES2.total.exp5.100.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp5.100.ICC.med[,j],RES2.total.exp5.100.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp5.100.ICC.med[,j],RES2.total.exp5.100.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp5.100.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp5.100.ICC.med<-rbind(result.data2.exp5.100.ICC.med,tot)
}


result.data2.exp5.100.ICC.med<-as.data.frame(result.data2.exp5.100.ICC.med)
colnames(result.data2.exp5.100.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.100.ICC.med$nselected<-as.numeric(as.character(result.data2.exp5.100.ICC.med$nselected))
result.data2.exp5.100.ICC.med$tn<-as.numeric(as.character(result.data2.exp5.100.ICC.med$tn))
result.data2.exp5.100.ICC.med$fn<-as.numeric(as.character(result.data2.exp5.100.ICC.med$fn))
result.data2.exp5.100.ICC.med$fp<-as.numeric(as.character(result.data2.exp5.100.ICC.med$fp))
result.data2.exp5.100.ICC.med$tp<-as.numeric(as.character(result.data2.exp5.100.ICC.med$tp))

result.data2.exp5.100.ICC.med$sensitivity<-round(result.data2.exp5.100.ICC.med$tp/(result.data2.exp5.100.ICC.med$tp+result.data2.exp5.100.ICC.med$fn)*100,1)
result.data2.exp5.100.ICC.med$specificity<-round(result.data2.exp5.100.ICC.med$tn/(result.data2.exp5.100.ICC.med$tn+result.data2.exp5.100.ICC.med$fp)*100,1)
result.data2.exp5.100.ICC.med$false.pos.rate<-round(result.data2.exp5.100.ICC.med$fp/(result.data2.exp5.100.ICC.med$fp+result.data2.exp5.100.ICC.med$tn)*100,1)
result.data2.exp5.100.ICC.med$false.neg.rate<-round(result.data2.exp5.100.ICC.med$fn/(result.data2.exp5.100.ICC.med$fn+result.data2.exp5.100.ICC.med$tp)*100,1)
result.data2.exp5.100.ICC.med$false.disc.rate<-round(result.data2.exp5.100.ICC.med$fp/(result.data2.exp5.100.ICC.med$fp+result.data2.exp5.100.ICC.med$tp)*100,1)
result.data2.exp5.100.ICC.med$false.omit.rate<-round(result.data2.exp5.100.ICC.med$fn/(result.data2.exp5.100.ICC.med$fn+result.data2.exp5.100.ICC.med$tn)*100,1)

result.data2.exp5.100.ICC.med$expo<-5
result.data2.exp5.100.ICC.med$ICC<-0.5


# High ICC
result.data2.exp5.100.ICC.high<-NULL
for (j in 4:(dim(RES2.total.exp5.100.ICC.high)[2]-1)){
  RES2.total.exp5.100.ICC.high[,j]<-as.factor(as.character(RES2.total.exp5.100.ICC.high[,j]))  
  RES2.total.exp5.100.ICC.high$true.pred<-as.factor(as.character(RES2.total.exp5.100.ICC.high$true.pred))
  tn<-squareTable(RES2.total.exp5.100.ICC.high[,j],RES2.total.exp5.100.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp5.100.ICC.high[,j],RES2.total.exp5.100.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp5.100.ICC.high[,j],RES2.total.exp5.100.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp5.100.ICC.high[,j],RES2.total.exp5.100.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp5.100.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp5.100.ICC.high<-rbind(result.data2.exp5.100.ICC.high,tot)
}


result.data2.exp5.100.ICC.high<-as.data.frame(result.data2.exp5.100.ICC.high)
colnames(result.data2.exp5.100.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp5.100.ICC.high$nselected<-as.numeric(as.character(result.data2.exp5.100.ICC.high$nselected))
result.data2.exp5.100.ICC.high$tn<-as.numeric(as.character(result.data2.exp5.100.ICC.high$tn))
result.data2.exp5.100.ICC.high$fn<-as.numeric(as.character(result.data2.exp5.100.ICC.high$fn))
result.data2.exp5.100.ICC.high$fp<-as.numeric(as.character(result.data2.exp5.100.ICC.high$fp))
result.data2.exp5.100.ICC.high$tp<-as.numeric(as.character(result.data2.exp5.100.ICC.high$tp))

result.data2.exp5.100.ICC.high$sensitivity<-round(result.data2.exp5.100.ICC.high$tp/(result.data2.exp5.100.ICC.high$tp+result.data2.exp5.100.ICC.high$fn)*100,1)
result.data2.exp5.100.ICC.high$specificity<-round(result.data2.exp5.100.ICC.high$tn/(result.data2.exp5.100.ICC.high$tn+result.data2.exp5.100.ICC.high$fp)*100,1)
result.data2.exp5.100.ICC.high$false.pos.rate<-round(result.data2.exp5.100.ICC.high$fp/(result.data2.exp5.100.ICC.high$fp+result.data2.exp5.100.ICC.high$tn)*100,1)
result.data2.exp5.100.ICC.high$false.neg.rate<-round(result.data2.exp5.100.ICC.high$fn/(result.data2.exp5.100.ICC.high$fn+result.data2.exp5.100.ICC.high$tp)*100,1)
result.data2.exp5.100.ICC.high$false.disc.rate<-round(result.data2.exp5.100.ICC.high$fp/(result.data2.exp5.100.ICC.high$fp+result.data2.exp5.100.ICC.high$tp)*100,1)
result.data2.exp5.100.ICC.high$false.omit.rate<-round(result.data2.exp5.100.ICC.high$fn/(result.data2.exp5.100.ICC.high$fn+result.data2.exp5.100.ICC.high$tn)*100,1)

result.data2.exp5.100.ICC.high$expo<-5
result.data2.exp5.100.ICC.high$ICC<-0.9


############### exp10

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/dataY2andX/exp10")
nexp<-10


load(file="RES2.all.exp10.100.2023.RData")
colnames(RES2.all.exp10.100)
colnames(RES2.all.exp10.100)<-c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                                "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")
RES2.all.exp10.100<-RES2.all.exp10.100[,c("numsim","var","true.pred","Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh","Raw.ExWAS.MLM.bh",
                                "Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA")]


load(file="twostep/RES2avRed.all.exp10.100.RData")
colnames(RES2avRed.all.exp10.100)
colnames(RES2avRed.all.exp10.100)<-c("numsim","var","true.pred","Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon","Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh",
                                     "Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA")


load(file="dlnm/RES2.DLNM.all.exp10.100.RData")
colnames(RES2.DLNM.all.exp10.100)
colnames(RES2.DLNM.all.exp10.100)<-c("numsim","var","true.pred","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                     "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                     "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",     
                                     "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")


load(file="dlnm/RES2.DLNM.AVG.all.exp10.100.RData")
colnames(RES2.DLNM.AVG.all.exp10.100)
colnames(RES2.DLNM.AVG.all.exp10.100)<-c("numsim","var","true.pred","Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by","Av.DLNMselect.none","Av.DLNMselect.bonf",
                                         "Av.DLNMselect.bh","Av.DLNMselect.by","Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by")


RES2.total.exp10.100<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c("var","numsim","true.pred")), 
                             list(RES2.all.exp10.100,RES2avRed.all.exp10.100,RES2.DLNM.all.exp10.100,RES2.DLNM.AVG.all.exp10.100))


### Merge avec ICC

ICC.data2.exp10.red<-matrix(ncol=dim(ICC.i[[1]])[2])
colnames(ICC.data2.exp10.red)<-colnames(ICC.i[[1]])
for (i in 1:nsim){
  ICC.data2.exp10.red<-rbind(ICC.data2.exp10.red,ICC.i[[i]])
}
ICC.data2.exp10.red<-ICC.data2.exp10.red[-1,]
ICC.data2.exp10.red$numsim<-as.integer(ICC.data2.exp10.red$numsim)
colnames(ICC.data2.exp10.red)<-c("var","numsim","ICC")

ICC.data2.exp10.red$expo_name<-ifelse(nchar(as.character(ICC.data2.exp10.red$var))==4,substring(ICC.data2.exp10.red$var,first=1,last=2),
                                      ifelse(nchar(as.character(ICC.data2.exp10.red$var))==5,substring(ICC.data2.exp10.red$var,first=1,last=3),
                                             substring(ICC.data2.exp10.red$var,first=1,last=4)))

ICC.data2.exp10.red<-ICC.data2.exp10.red[,c("expo_name","numsim","ICC")]
ICC.data2.exp10.red<-unique(ICC.data2.exp10.red)
colnames(ICC.data2.exp10.red)<-c("var","numsim","ICC")

RES2.total.exp10.100.ICC<-merge(RES2.total.exp10.100,ICC.data2.exp10.red,by=c("var","numsim"))

table(RES2.total.exp10.100.ICC$true.pred,RES2.total.exp10.100.ICC$ICC)


squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


RES2.total.exp10.100.ICC.low<-RES2.total.exp10.100.ICC[RES2.total.exp10.100.ICC$ICC==0.1,]
RES2.total.exp10.100.ICC.med<-RES2.total.exp10.100.ICC[RES2.total.exp10.100.ICC$ICC==0.5,]
RES2.total.exp10.100.ICC.high<-RES2.total.exp10.100.ICC[RES2.total.exp10.100.ICC$ICC==0.9,]


### All
RES2.i.total.exp10.100.ICC<-vector("list", nsim)
for (i in 1:nsim){
  RES2.i.total.exp10.100.ICC[[i]]<-RES2.total.exp10.100.ICC[RES2.total.exp10.100.ICC$numsim == i,]
}

result.data2.i.exp10.100<- vector("list", nsim)
for (i in 1:nsim){
  for (j in 4:(dim(RES2.i.total.exp10.100.ICC[[i]])[2]-1)){
    RES2.i.total.exp10.100.ICC[[i]][,j]<-as.factor(as.character(RES2.i.total.exp10.100.ICC[[i]][,j]))  
    RES2.i.total.exp10.100.ICC[[i]]$true.pred<-as.factor(as.character(RES2.i.total.exp10.100.ICC[[i]]$true.pred))
    tn<-squareTable(RES2.i.total.exp10.100.ICC[[i]][,j],RES2.i.total.exp10.100.ICC[[i]]$true.pred)[1,1]
    fn<-squareTable(RES2.i.total.exp10.100.ICC[[i]][,j],RES2.i.total.exp10.100.ICC[[i]]$true.pred)[1,2]
    fp<-squareTable(RES2.i.total.exp10.100.ICC[[i]][,j],RES2.i.total.exp10.100.ICC[[i]]$true.pred)[2,1]
    tp<-squareTable(RES2.i.total.exp10.100.ICC[[i]][,j],RES2.i.total.exp10.100.ICC[[i]]$true.pred)[2,2]
    nselected<-tp+fp
    tot<-c(colnames(RES2.i.total.exp10.100.ICC[[i]][j]),nselected,tn,fn,fp,tp,i)
    result.data2.i.exp10.100[[i]]<-rbind(result.data2.i.exp10.100[[i]],tot)
  }
}


result.data2.exp10.100<-matrix(ncol=dim(result.data2.i.exp10.100[[1]])[2])
for (i in 1:nsim){
  result.data2.exp10.100<-rbind(result.data2.exp10.100,result.data2.i.exp10.100[[i]])
}

result.data2.exp10.100<-as.data.frame(result.data2.exp10.100)
colnames(result.data2.exp10.100)<-c("Method","nselected","tn","fn","fp","tp","data.i")
result.data2.exp10.100$nselected<-as.numeric(as.character(result.data2.exp10.100$nselected))
result.data2.exp10.100$tn<-as.numeric(as.character(result.data2.exp10.100$tn))
result.data2.exp10.100$fn<-as.numeric(as.character(result.data2.exp10.100$fn))
result.data2.exp10.100$fp<-as.numeric(as.character(result.data2.exp10.100$fp))
result.data2.exp10.100$tp<-as.numeric(as.character(result.data2.exp10.100$tp))

result.data2.exp10.100$sensitivity<-round(result.data2.exp10.100$tp/(result.data2.exp10.100$tp+result.data2.exp10.100$fn)*100,1)
result.data2.exp10.100$specificity<-round(result.data2.exp10.100$tn/(result.data2.exp10.100$tn+result.data2.exp10.100$fp)*100,1)
result.data2.exp10.100$false.pos.rate<-round(result.data2.exp10.100$fp/(result.data2.exp10.100$fp+result.data2.exp10.100$tn)*100,1)
result.data2.exp10.100$false.neg.rate<-round(result.data2.exp10.100$fn/(result.data2.exp10.100$fn+result.data2.exp10.100$tp)*100,1)
result.data2.exp10.100$false.disc.rate<-round(result.data2.exp10.100$fp/(result.data2.exp10.100$fp+result.data2.exp10.100$tp)*100,1)
result.data2.exp10.100$false.omit.rate<-round(result.data2.exp10.100$fn/(result.data2.exp10.100$fn+result.data2.exp10.100$tn)*100,1)

result.data2.exp10.100$expo<-10
result.data2.exp10.100$ICC<-1


# Low ICC
result.data2.exp10.100.ICC.low<-NULL
for (j in 4:(dim(RES2.total.exp10.100.ICC.low)[2]-1)){
  RES2.total.exp10.100.ICC.low[,j]<-as.factor(as.character(RES2.total.exp10.100.ICC.low[,j]))  
  RES2.total.exp10.100.ICC.low$true.pred<-as.factor(as.character(RES2.total.exp10.100.ICC.low$true.pred))
  tn<-squareTable(RES2.total.exp10.100.ICC.low[,j],RES2.total.exp10.100.ICC.low$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp10.100.ICC.low[,j],RES2.total.exp10.100.ICC.low$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp10.100.ICC.low[,j],RES2.total.exp10.100.ICC.low$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp10.100.ICC.low[,j],RES2.total.exp10.100.ICC.low$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp10.100.ICC.low[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp10.100.ICC.low<-rbind(result.data2.exp10.100.ICC.low,tot)
}


result.data2.exp10.100.ICC.low<-as.data.frame(result.data2.exp10.100.ICC.low)
colnames(result.data2.exp10.100.ICC.low)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10.100.ICC.low$nselected<-as.numeric(as.character(result.data2.exp10.100.ICC.low$nselected))
result.data2.exp10.100.ICC.low$tn<-as.numeric(as.character(result.data2.exp10.100.ICC.low$tn))
result.data2.exp10.100.ICC.low$fn<-as.numeric(as.character(result.data2.exp10.100.ICC.low$fn))
result.data2.exp10.100.ICC.low$fp<-as.numeric(as.character(result.data2.exp10.100.ICC.low$fp))
result.data2.exp10.100.ICC.low$tp<-as.numeric(as.character(result.data2.exp10.100.ICC.low$tp))

result.data2.exp10.100.ICC.low$sensitivity<-round(result.data2.exp10.100.ICC.low$tp/(result.data2.exp10.100.ICC.low$tp+result.data2.exp10.100.ICC.low$fn)*100,1)
result.data2.exp10.100.ICC.low$specificity<-round(result.data2.exp10.100.ICC.low$tn/(result.data2.exp10.100.ICC.low$tn+result.data2.exp10.100.ICC.low$fp)*100,1)
result.data2.exp10.100.ICC.low$false.pos.rate<-round(result.data2.exp10.100.ICC.low$fp/(result.data2.exp10.100.ICC.low$fp+result.data2.exp10.100.ICC.low$tn)*100,1)
result.data2.exp10.100.ICC.low$false.neg.rate<-round(result.data2.exp10.100.ICC.low$fn/(result.data2.exp10.100.ICC.low$fn+result.data2.exp10.100.ICC.low$tp)*100,1)
result.data2.exp10.100.ICC.low$false.disc.rate<-round(result.data2.exp10.100.ICC.low$fp/(result.data2.exp10.100.ICC.low$fp+result.data2.exp10.100.ICC.low$tp)*100,1)
result.data2.exp10.100.ICC.low$false.omit.rate<-round(result.data2.exp10.100.ICC.low$fn/(result.data2.exp10.100.ICC.low$fn+result.data2.exp10.100.ICC.low$tn)*100,1)

result.data2.exp10.100.ICC.low$expo<-10
result.data2.exp10.100.ICC.low$ICC<-0.1



# Medium ICC
result.data2.exp10.100.ICC.med<-NULL
for (j in 4:(dim(RES2.total.exp10.100.ICC.med)[2]-1)){
  RES2.total.exp10.100.ICC.med[,j]<-as.factor(as.character(RES2.total.exp10.100.ICC.med[,j]))  
  RES2.total.exp10.100.ICC.med$true.pred<-as.factor(as.character(RES2.total.exp10.100.ICC.med$true.pred))
  tn<-squareTable(RES2.total.exp10.100.ICC.med[,j],RES2.total.exp10.100.ICC.med$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp10.100.ICC.med[,j],RES2.total.exp10.100.ICC.med$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp10.100.ICC.med[,j],RES2.total.exp10.100.ICC.med$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp10.100.ICC.med[,j],RES2.total.exp10.100.ICC.med$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp10.100.ICC.med[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp10.100.ICC.med<-rbind(result.data2.exp10.100.ICC.med,tot)
}


result.data2.exp10.100.ICC.med<-as.data.frame(result.data2.exp10.100.ICC.med)
colnames(result.data2.exp10.100.ICC.med)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10.100.ICC.med$nselected<-as.numeric(as.character(result.data2.exp10.100.ICC.med$nselected))
result.data2.exp10.100.ICC.med$tn<-as.numeric(as.character(result.data2.exp10.100.ICC.med$tn))
result.data2.exp10.100.ICC.med$fn<-as.numeric(as.character(result.data2.exp10.100.ICC.med$fn))
result.data2.exp10.100.ICC.med$fp<-as.numeric(as.character(result.data2.exp10.100.ICC.med$fp))
result.data2.exp10.100.ICC.med$tp<-as.numeric(as.character(result.data2.exp10.100.ICC.med$tp))

result.data2.exp10.100.ICC.med$sensitivity<-round(result.data2.exp10.100.ICC.med$tp/(result.data2.exp10.100.ICC.med$tp+result.data2.exp10.100.ICC.med$fn)*100,1)
result.data2.exp10.100.ICC.med$specificity<-round(result.data2.exp10.100.ICC.med$tn/(result.data2.exp10.100.ICC.med$tn+result.data2.exp10.100.ICC.med$fp)*100,1)
result.data2.exp10.100.ICC.med$false.pos.rate<-round(result.data2.exp10.100.ICC.med$fp/(result.data2.exp10.100.ICC.med$fp+result.data2.exp10.100.ICC.med$tn)*100,1)
result.data2.exp10.100.ICC.med$false.neg.rate<-round(result.data2.exp10.100.ICC.med$fn/(result.data2.exp10.100.ICC.med$fn+result.data2.exp10.100.ICC.med$tp)*100,1)
result.data2.exp10.100.ICC.med$false.disc.rate<-round(result.data2.exp10.100.ICC.med$fp/(result.data2.exp10.100.ICC.med$fp+result.data2.exp10.100.ICC.med$tp)*100,1)
result.data2.exp10.100.ICC.med$false.omit.rate<-round(result.data2.exp10.100.ICC.med$fn/(result.data2.exp10.100.ICC.med$fn+result.data2.exp10.100.ICC.med$tn)*100,1)

result.data2.exp10.100.ICC.med$expo<-10
result.data2.exp10.100.ICC.med$ICC<-0.5


# High ICC
result.data2.exp10.100.ICC.high<-NULL
for (j in 4:(dim(RES2.total.exp10.100.ICC.high)[2]-1)){
  RES2.total.exp10.100.ICC.high[,j]<-as.factor(as.character(RES2.total.exp10.100.ICC.high[,j]))  
  RES2.total.exp10.100.ICC.high$true.pred<-as.factor(as.character(RES2.total.exp10.100.ICC.high$true.pred))
  tn<-squareTable(RES2.total.exp10.100.ICC.high[,j],RES2.total.exp10.100.ICC.high$true.pred)[1,1]
  fn<-squareTable(RES2.total.exp10.100.ICC.high[,j],RES2.total.exp10.100.ICC.high$true.pred)[1,2]
  fp<-squareTable(RES2.total.exp10.100.ICC.high[,j],RES2.total.exp10.100.ICC.high$true.pred)[2,1]
  tp<-squareTable(RES2.total.exp10.100.ICC.high[,j],RES2.total.exp10.100.ICC.high$true.pred)[2,2]
  nselected<-tp+fp
  tot<-c(colnames(RES2.total.exp10.100.ICC.high[j]),nselected,tn,fn,fp,tp,i)
  result.data2.exp10.100.ICC.high<-rbind(result.data2.exp10.100.ICC.high,tot)
}


result.data2.exp10.100.ICC.high<-as.data.frame(result.data2.exp10.100.ICC.high)
colnames(result.data2.exp10.100.ICC.high)<-c("Method","nselected","tn","fn","fp","tp","data.i")

result.data2.exp10.100.ICC.high$nselected<-as.numeric(as.character(result.data2.exp10.100.ICC.high$nselected))
result.data2.exp10.100.ICC.high$tn<-as.numeric(as.character(result.data2.exp10.100.ICC.high$tn))
result.data2.exp10.100.ICC.high$fn<-as.numeric(as.character(result.data2.exp10.100.ICC.high$fn))
result.data2.exp10.100.ICC.high$fp<-as.numeric(as.character(result.data2.exp10.100.ICC.high$fp))
result.data2.exp10.100.ICC.high$tp<-as.numeric(as.character(result.data2.exp10.100.ICC.high$tp))

result.data2.exp10.100.ICC.high$sensitivity<-round(result.data2.exp10.100.ICC.high$tp/(result.data2.exp10.100.ICC.high$tp+result.data2.exp10.100.ICC.high$fn)*100,1)
result.data2.exp10.100.ICC.high$specificity<-round(result.data2.exp10.100.ICC.high$tn/(result.data2.exp10.100.ICC.high$tn+result.data2.exp10.100.ICC.high$fp)*100,1)
result.data2.exp10.100.ICC.high$false.pos.rate<-round(result.data2.exp10.100.ICC.high$fp/(result.data2.exp10.100.ICC.high$fp+result.data2.exp10.100.ICC.high$tn)*100,1)
result.data2.exp10.100.ICC.high$false.neg.rate<-round(result.data2.exp10.100.ICC.high$fn/(result.data2.exp10.100.ICC.high$fn+result.data2.exp10.100.ICC.high$tp)*100,1)
result.data2.exp10.100.ICC.high$false.disc.rate<-round(result.data2.exp10.100.ICC.high$fp/(result.data2.exp10.100.ICC.high$fp+result.data2.exp10.100.ICC.high$tp)*100,1)
result.data2.exp10.100.ICC.high$false.omit.rate<-round(result.data2.exp10.100.ICC.high$fn/(result.data2.exp10.100.ICC.high$fn+result.data2.exp10.100.ICC.high$tn)*100,1)

result.data2.exp10.100.ICC.high$expo<-10
result.data2.exp10.100.ICC.high$ICC<-0.9


result.data2.100.ICC<-rbind(result.data2.exp3.100,result.data2.exp3.100.ICC.low,result.data2.exp3.100.ICC.med,result.data2.exp3.100.ICC.high,
                        result.data2.exp5.100,result.data2.exp5.100.ICC.low,result.data2.exp5.100.ICC.med,result.data2.exp5.100.ICC.high,
                        result.data2.exp10.100,result.data2.exp10.100.ICC.low,result.data2.exp10.100.ICC.med,result.data2.exp10.100.ICC.high)

result.data2.100.ICC$false.disc.rate[result.data2.100.ICC$nselected==0]<-0

save(result.data2.100.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data2.100_ICC.2023.Rdata")
write.csv2(result.data2.100.ICC,file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary.100_data2_ICC.2023.csv")


