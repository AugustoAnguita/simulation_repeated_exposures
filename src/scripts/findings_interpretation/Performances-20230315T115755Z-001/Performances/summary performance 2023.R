library(ggplot2)
library(dplyr)
library(corrplot)

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/")



############################################## correlation matrix between baseline exposures ###########################################


load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/baseline/resu.sim.dataX.i.RData")

# 1st imputed dataset // all variables n=500
expo1<-resu.sim.dataX.i[[1]]$X
cor_mat<-cor(expo1)
corrplot(cor_mat, tl.pos='n')
corrplot(cor_mat, order = "hclust", tl.pos='n')

# 1st imputed dataset // first time point n=100
expo.red<-resu.sim.dataX.i[[1]]$X
retain<-seq(1,500,by=5)
expo.red<-expo.red[,retain]
cor_mat<-cor(expo.red)
corrplot(cor_mat, tl.pos='n')
corrplot(cor_mat, order = "hclust", tl.pos='n')


# 1st imputed dataset // first variables n=100
expo.red<-resu.sim.dataX.i[[1]]$X
cor_mat<-cor(expo.red[1:100])
corrplot(cor_mat, tl.pos='n')
corrplot(cor_mat, order = "hclust", tl.pos='n')


result.perf.expo.dlnm.exp5.100




################################### Performance to identify the true exposure whatever the time point ####################################



######################################
###             DATA1              ###
######################################


# data1 = all time points - raw data

load(file="dataY1andX/exp3/result.perf.expo.exp3.2023.RData")
result.perf.expo.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/result.perf.expo.exp5.2023.RData")
result.perf.expo.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/result.perf.expo.exp10.2023.RData")
result.perf.expo.exp10$nb.pred<-"10"

result.data1.100<-rbind(result.perf.expo.exp3,result.perf.expo.exp5,result.perf.expo.exp10)
colnames(result.data1.100)

result.data1.100<-result.data1.100[result.data1.100$Method %in% c("EWAS.TP.none_sum", "EWAS.TP.bon_sum", "EWAS.TP.bh_sum", "EWAS.TP.by_sum", 
                                                                  "EWAS_LM.TP.none_sum", "EWAS_LM.TP.bon_sum", "EWAS_LM.TP.bh_sum", "EWAS_LM.TP.by_sum",
                                                                  "ENET_MIN.TP_sum","ENET_OPT.TP_sum","SPLS_MIN.TP_sum","sNPLS.TP_sum","MMPC.TP_sum","DSA.TP_sum"),]

result.data1.100$Method <- factor(result.data1.100$Method, levels = c("EWAS.TP.none_sum", "EWAS.TP.bon_sum", "EWAS.TP.bh_sum", "EWAS.TP.by_sum", 
                                                                      "EWAS_LM.TP.none_sum", "EWAS_LM.TP.bon_sum", "EWAS_LM.TP.bh_sum", "EWAS_LM.TP.by_sum",
                                                                      "ENET_MIN.TP_sum","ENET_OPT.TP_sum","SPLS_MIN.TP_sum","sNPLS.TP_sum","MMPC.TP_sum","DSA.TP_sum"))

levels(result.data1.100$Method) <- list(Raw.EWAS.none="EWAS.TP.none_sum", Raw.EWAS.bon="EWAS.TP.bon_sum", Raw.EWAS.bh="EWAS.TP.bh_sum", Raw.EWAS.by="EWAS.TP.by_sum", 
                                        Raw.EWAS_LM.none="EWAS_LM.TP.none_sum", Raw.EWAS_LM.bon="EWAS_LM.TP.bon_sum", Raw.EWAS_LM.bh="EWAS_LM.TP.bh_sum", Raw.EWAS_LM.by="EWAS_LM.TP.by_sum",
                                        Raw.ENET_MIN="ENET_MIN.TP_sum",Raw.ENET_OPT="ENET_OPT.TP_sum",Raw.SPLS="SPLS_MIN.TP_sum",Raw.sNPLS="sNPLS.TP_sum",Raw.MMPC="MMPC.TP_sum",Raw.DSA="DSA.TP_sum")

result.data1.100$nb.pred <- factor(result.data1.100$nb.pred, levels = c("3","5","10"))


# data1 = all time points - average - results obtained after the 1st step (=mean_variable)

#### the following datasets corresponds to the performance to identify the true expo at the 1-step (not after the 2-step) approach applied with all methods
load(file="dataY1andX/exp3/twostep/result.data1av.exp3.100.RData")
result.data1av.exp3.100$nb.pred<-"3"
load(file="dataY1andX/exp5/twostep/result.data1av.exp5.100.RData")
result.data1av.exp5.100$nb.pred<-"5"
load(file="dataY1andX/exp10/twostep/result.data1av.exp10.100.RData")
result.data1av.exp10.100$nb.pred<-"10"

result.data1av.100<-rbind(result.data1av.exp3.100,result.data1av.exp5.100,result.data1av.exp10.100)


colnames(result.data1av.100)
unique(result.data1av.100$Method)

result.data1av.100$Method <- factor(result.data1av.100$Method, levels = c("EWAS.TP.none", "EWAS.TP.bon", "EWAS.TP.bh", 
                                                                          "EWAS.TP.by","EWAS_LM.TP.none", "EWAS_LM.TP.bon",
                                                                          "EWAS_LM.TP.bh", "EWAS_LM.TP.by", "ENET_MIN.TP",
                                                                          "ENET_OPT.TP","SPLS_MIN.TP","MMPC.TP","DSA.TP"))

levels(result.data1av.100$Method) <- list(Av.EWAS.none="EWAS.TP.none", Av.EWAS.bon="EWAS.TP.bon", Av.EWAS.bh="EWAS.TP.bh", Av.EWAS.by="EWAS.TP.by", 
                                          Av.EWAS_LM.none="EWAS_LM.TP.none", Av.EWAS_LM.bon="EWAS_LM.TP.bon", Av.EWAS_LM.bh="EWAS_LM.TP.bh", Av.EWAS_LM.by="EWAS_LM.TP.by",
                                          Av.ENET_MIN="ENET_MIN.TP",Av.ENET_OPT="ENET_OPT.TP",Av.SPLS="SPLS_MIN.TP",Av.MMPC="MMPC.TP",Av.DSA="DSA.TP")


result.data1av.100$nb.pred <- factor(result.data1av.100$nb.pred, levels = c("3","5","10"))



# data1 = all time points - average - results obtained after the 2nd step (=raw_variable)

#### the following datasets corresponds to the performance to identify the true expo at the 2-step (=any raw variable significant within group of 5)
load(file="dataY1andX/exp3/twostep/result.perf.expo.av.exp3.RData")
result.perf.expo.av.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/twostep/result.perf.expo.av.exp5.RData")
result.perf.expo.av.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/twostep/result.perf.expo.av.exp10.RData")
result.perf.expo.av.exp10$nb.pred<-"10"

result.data1av.2step.100<-rbind(result.perf.expo.av.exp3,result.perf.expo.av.exp5,result.perf.expo.av.exp10)


colnames(result.data1av.2step.100)
unique(result.data1av.2step.100$Method)

result.data1av.2step.100$Method <- factor(result.data1av.2step.100$Method, levels = c("EWAS.TP.none_sum", "EWAS.TP.bon_sum", "EWAS.TP.bh_sum", 
                                                                                      "EWAS.TP.by_sum","EWAS_LM.TP.none_sum", "EWAS_LM.TP.bon_sum",
                                                                                      "EWAS_LM.TP.bh_sum", "EWAS_LM.TP.by_sum", "ENET_MIN.TP_sum",
                                                                                      "ENET_OPT.TP_sum","SPLS_MIN.TP_sum","MMPC.TP_sum","DSA.TP_sum"))

levels(result.data1av.2step.100$Method) <- list(Av.2step.EWAS.none="EWAS.TP.none_sum", Av.2step.EWAS.bon="EWAS.TP.bon_sum", Av.2step.EWAS.bh="EWAS.TP.bh_sum", Av.2step.EWAS.by="EWAS.TP.by_sum", 
                                                Av.2step.EWAS_LM.none="EWAS_LM.TP.none_sum", Av.2step.EWAS_LM.bon="EWAS_LM.TP.bon_sum", Av.2step.EWAS_LM.bh="EWAS_LM.TP.bh_sum", Av.2step.EWAS_LM.by="EWAS_LM.TP.by_sum",
                                                Av.2step.ENET_MIN="ENET_MIN.TP_sum",Av.2step.ENET_OPT="ENET_OPT.TP_sum",Av.2step.SPLS="SPLS_MIN.TP_sum",Av.2step.MMPC="MMPC.TP_sum",Av.2step.DSA="DSA.TP_sum")


result.data1av.2step.100$nb.pred <- factor(result.data1av.2step.100$nb.pred, levels = c("3","5","10"))




# data1 - all time point - DLNM - 1st step (crossbasis)

load(file="dataY1andX/exp3/dlnm/result.perf.expo.dlnm.exp3.100.RData")
result.perf.expo.dlnm.exp3.100$nb.pred<-"3"
load(file="dataY1andX/exp5/dlnm/result.perf.expo.dlnm.exp5.100.RData")
result.perf.expo.dlnm.exp5.100$nb.pred<-"5"
load(file="dataY1andX/exp10/dlnm/result.perf.expo.dlnm.exp10.100.RData")
result.perf.expo.dlnm.exp10.100$nb.pred<-"10"

result.data1.dlnm.100<-rbind(result.perf.expo.dlnm.exp3.100,result.perf.expo.dlnm.exp5.100,result.perf.expo.dlnm.exp10.100)

colnames(result.data1.dlnm.100)
unique(result.data1.dlnm.100$Method)

result.data1.dlnm.100$Method <- factor(result.data1.dlnm.100$Method, levels = c("DLNMpen.none.TP01.100_sum","DLNMpen.bonf.TP01.100_sum","DLNMpen.bh.TP01.100_sum","DLNMpen.by.TP01.100_sum",
                                                                                "DLNMselect.none.TP01.100_sum","DLNMselect.bonf.TP01.100_sum","DLNMselect.bh.TP01.100_sum","DLNMselect.by.TP01.100_sum",
                                                                                "DLNMselectback.none.TP01.100_sum","DLNMselectback.bonf.TP01.100_sum","DLNMselectback.bh.TP01.100_sum","DLNMselectback.by.TP01.100_sum",
                                                                                "DLNMselectforward.none.TP01.100_sum","DLNMselectforward.bonf.TP01.100_sum","DLNMselectforward.bh.TP01.100_sum","DLNMselectforward.by.TP01.100_sum"))

levels(result.data1.dlnm.100$Method) <- list(DLNMpen.none="DLNMpen.none.TP01.100_sum", DLNMpen.bonf="DLNMpen.bonf.TP01.100_sum", DLNMpen.bh="DLNMpen.bh.TP01.100_sum", DLNMpen.by="DLNMpen.by.TP01.100_sum", 
                                             DLNMselect.none="DLNMselect.none.TP01.100_sum", DLNMselect.bonf="DLNMselect.bonf.TP01.100_sum", DLNMselect.bh="DLNMselect.bh.TP01.100_sum", DLNMselect.by="DLNMselect.by.TP01.100_sum",
                                             DLNMselectback.none="DLNMselectback.none.TP01.100_sum",DLNMselectback.bonf="DLNMselectback.bonf.TP01.100_sum",DLNMselectback.bh="DLNMselectback.bh.TP01.100_sum",DLNMselectback.by="DLNMselectback.by.TP01.100_sum",
                                             DLNMselectforward.none="DLNMselectforward.none.TP01.100_sum",DLNMselectforward.bonf="DLNMselectforward.bonf.TP01.100_sum",DLNMselectforward.bh="DLNMselectforward.bh.TP01.100_sum",DLNMselectforward.by="DLNMselectforward.by.TP01.100_sum")

result.data1.dlnm.100$nb.pred <- factor(result.data1.dlnm.100$nb.pred, levels = c("3","5","10"))

result.data1.dlnm.100<-result.data1.dlnm.100[! result.data1.dlnm.100$Method %in% c("DLNMselectforward.bonf", "DLNMselectforward.bh", "DLNMselectforward.by"),]

# data1 - all time point - DLNM - 2nd step (lag)

load(file="dataY1andX/exp3/dlnm/result.perf.expo.dlnm.exp3.RData")
result.perf.expo.dlnm.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/dlnm/result.perf.expo.dlnm.exp5.RData")
result.perf.expo.dlnm.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/dlnm/result.perf.expo.dlnm.exp10.RData")
result.perf.expo.dlnm.exp10$nb.pred<-"10"

result.data1.dlnm.step2.100<-rbind(result.perf.expo.dlnm.exp3,result.perf.expo.dlnm.exp5,result.perf.expo.dlnm.exp10)

colnames(result.data1.dlnm.step2.100)
unique(result.data1.dlnm.step2.100$Method)

result.data1.dlnm.step2.100$Method <- factor(result.data1.dlnm.step2.100$Method, levels = c("DLNMpen.none_sum","DLNMpen.bonf_sum","DLNMpen.bh_sum","DLNMpen.by_sum",
                                                                                            "DLNMselect.none_sum","DLNMselect.bonf_sum","DLNMselect.bh_sum","DLNMselect.by_sum",
                                                                                            "DLNMselectback.none_sum","DLNMselectback.bonf_sum","DLNMselectback.bh_sum","DLNMselectback.by_sum",
                                                                                            "DLNMselectforward.none_sum","DLNMselectforward.bonf_sum","DLNMselectforward.bh_sum","DLNMselectforward.by_sum"))

levels(result.data1.dlnm.step2.100$Method) <- list(DLNMpen.2step.none="DLNMpen.none_sum", DLNMpen.2step.bonf="DLNMpen.bonf_sum", DLNMpen.2step.bh="DLNMpen.bh_sum", DLNMpen.2step.by="DLNMpen.by_sum", 
                                                   DLNMselect.2step.none="DLNMselect.none_sum", DLNMselect.2step.bonf="DLNMselect.bonf_sum", DLNMselect.2step.bh="DLNMselect.bh_sum", DLNMselect.2step.by="DLNMselect.by_sum",
                                                   DLNMselectback.2step.none="DLNMselectback.none_sum",DLNMselectback.2step.bonf="DLNMselectback.bonf_sum",DLNMselectback.2step.bh="DLNMselectback.bh_sum",DLNMselectback.2step.by="DLNMselectback.by_sum",
                                                   DLNMselectforward.2step.none="DLNMselectforward.none_sum",DLNMselectforward.2step.bonf="DLNMselectforward.bonf_sum",DLNMselectforward.2step.bh="DLNMselectforward.bh_sum",DLNMselectforward.2step.by="DLNMselectforward.by_sum")

result.data1.dlnm.step2.100$nb.pred <- factor(result.data1.dlnm.step2.100$nb.pred, levels = c("3","5","10"))

result.data1.dlnm.step2.100<-result.data1.dlnm.step2.100[! result.data1.dlnm.step2.100$Method %in% c("DLNMselectforward.2step.bonf", "DLNMselectforward.2step.bh", "DLNMselectforward.2step.by"),]


# data1 - all time point - DLNM AVG - 1st step (crossbasis)

load(file="dataY1andX/exp3/dlnm/result.perf.expo.dlnm.AVG.exp3.100.RData")
result.perf.expo.dlnm.AVG.exp3.100$nb.pred<-"3"
load(file="dataY1andX/exp5/dlnm/result.perf.expo.dlnm.AVG.exp5.100.RData")
result.perf.expo.dlnm.AVG.exp5.100$nb.pred<-"5"
load(file="dataY1andX/exp10/dlnm/result.perf.expo.dlnm.AVG.exp10.100.RData")
result.perf.expo.dlnm.AVG.exp10.100$nb.pred<-"10"

result.data1.dlnm.AVG.100<-rbind(result.perf.expo.dlnm.AVG.exp3.100,result.perf.expo.dlnm.AVG.exp5.100,result.perf.expo.dlnm.AVG.exp10.100)

colnames(result.data1.dlnm.AVG.100)
unique(result.data1.dlnm.AVG.100$Method)

result.data1.dlnm.AVG.100$Method <- factor(result.data1.dlnm.AVG.100$Method, levels = c("DLNMpen.TP01.100_sum","DLNMselect.TP01.100_sum","DLNMselectback.TP01.100_sum"))

levels(result.data1.dlnm.AVG.100$Method) <- list(DLNMpen.AVG="DLNMpen.TP01.100_sum", DLNMselect.AVG="DLNMselect.TP01.100_sum", DLNMselectback.AVG="DLNMselectback.TP01.100_sum")

result.data1.dlnm.AVG.100$nb.pred <- factor(result.data1.dlnm.AVG.100$nb.pred, levels = c("3","5","10"))


# data1 - all time point - DLNM AVG - 2nd step (lag)

load(file="dataY1andX/exp3/dlnm/result.perf.expo.dlnm.AVG.exp3.RData")
result.perf.expo.dlnm.AVG.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/dlnm/result.perf.expo.dlnm.AVG.exp5.RData")
result.perf.expo.dlnm.AVG.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/dlnm/result.perf.expo.dlnm.AVG.exp10.RData")
result.perf.expo.dlnm.AVG.exp10$nb.pred<-"10"

result.data1.dlnm.AVG.step2.100<-rbind(result.perf.expo.dlnm.AVG.exp3,result.perf.expo.dlnm.AVG.exp5,result.perf.expo.dlnm.AVG.exp10)

colnames(result.data1.dlnm.AVG.step2.100)
unique(result.data1.dlnm.AVG.step2.100$Method)

result.data1.dlnm.AVG.step2.100$Method <- factor(result.data1.dlnm.AVG.step2.100$Method, levels = c("DLNMpen.none_sum","DLNMpen.bonf_sum","DLNMpen.bh_sum","DLNMpen.by_sum",
                                                                                                    "DLNMselect.none_sum","DLNMselect.bonf_sum","DLNMselect.bh_sum","DLNMselect.by_sum",
                                                                                                    "DLNMselectback.none_sum","DLNMselectback.bonf_sum","DLNMselectback.bh_sum","DLNMselectback.by_sum"))

levels(result.data1.dlnm.AVG.step2.100$Method) <- list(DLNMpen.AVG.2step.none="DLNMpen.none_sum", DLNMpen.AVG.2step.bonf="DLNMpen.bonf_sum", DLNMpen.AVG.2step.bh="DLNMpen.bh_sum", DLNMpen.AVG.2step.by="DLNMpen.by_sum", 
                                                       DLNMselect.AVG.2step.none="DLNMselect.none_sum", DLNMselect.AVG.2step.bonf="DLNMselect.bonf_sum", DLNMselect.AVG.2step.bh="DLNMselect.bh_sum", DLNMselect.AVG.2step.by="DLNMselect.by_sum",
                                                       DLNMselectback.AVG.2step.none="DLNMselectback.none_sum",DLNMselectback.AVG.2step.bonf="DLNMselectback.bonf_sum",DLNMselectback.AVG.2step.bh="DLNMselectback.bh_sum",DLNMselectback.AVG.2step.by="DLNMselectback.by_sum")

result.data1.dlnm.AVG.step2.100$nb.pred <- factor(result.data1.dlnm.AVG.step2.100$nb.pred, levels = c("3","5","10"))




# total
data1.all.100<-rbind(result.data1.100,result.data1av.100,result.data1av.2step.100,result.data1.dlnm.100,result.data1.dlnm.step2.100,result.data1.dlnm.AVG.100,result.data1.dlnm.AVG.step2.100)

save(data1.all.100,file="Summary/data1.all.100.2023.Rdata")

data1.100<-data1.all.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity),
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate)
  )


data1.100$Sensitivity<-paste0(round(data1.100$sensitivity_mean,1)," (",round(data1.100$sensitivity_sd,1),")")
data1.100$FDR<-paste0(round(data1.100$false.disc.rate_mean,1)," (",round(data1.100$false.disc.rate_sd,1),")")


library(reshape2)
summary_data1.100<-dcast(melt(data1.100, id.vars=c("Method", "nb.pred")), Method~variable+nb.pred)

write.csv2(summary_data1.100,file="Summary/summary.data1.100.2023.csv")
save(summary_data1.100,file="Summary/summary.data1.100.2023.Rdata")



#################### plot


data<-result.data1.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data$Method)


data$color.perf<-as.factor(ifelse(data$sensitivity_mean>80 & data$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                  ifelse(data$sensitivity_mean>70 & data$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                         ifelse(data$sensitivity_mean>60 & data$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))

# raw data
plot1.raw<-data[data$Method%in%c("Raw.EWAS.none","Raw.EWAS.bon","Raw.EWAS.bh","Raw.EWAS.by","Raw.EWAS_LM.none","Raw.EWAS_LM.bon",
                                 "Raw.EWAS_LM.bh","Raw.EWAS_LM.by","Raw.ENET_MIN","Raw.ENET_OPT","Raw.SPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

ggplot(data = plot1.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1: All time points associated with Y - Raw data - Performance to identify the true exposure whatever the true time point")


# averaged levels - 1st step

data.av<-result.data1av.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.av$Method)

data.av$color.perf<-as.factor(ifelse(data.av$sensitivity_mean>80 & data.av$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                     ifelse(data.av$sensitivity_mean>70 & data.av$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                            ifelse(data.av$sensitivity_mean>60 & data.av$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot1.av<-data.av[data.av$Method%in%c("Av.EWAS.none","Av.EWAS.bon","Av.EWAS.bh","Av.EWAS.by","Av.EWAS_LM.none","Av.EWAS_LM.bon",
                                      "Av.EWAS_LM.bh","Av.EWAS_LM.by","Av.ENET_MIN","Av.ENET_OPT","Av.SPLS","Av.MMPC","Av.DSA"),]


ggplot(data = plot1.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1 - Performance to identify the true exposure whatever the true time point - Averaged exposure levels")




# averaged levels - 2-step

data.av.2step<-result.data1av.2step.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.av.2step$Method)

data.av.2step$color.perf<-as.factor(ifelse(data.av.2step$sensitivity_mean>80 & data.av.2step$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                           ifelse(data.av.2step$sensitivity_mean>70 & data.av.2step$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                                  ifelse(data.av.2step$sensitivity_mean>60 & data.av.2step$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot1.av.2step<-data.av.2step[data.av.2step$Method%in%c("Av.2step.EWAS.none","Av.2step.EWAS.bon","Av.2step.EWAS.bh","Av.2step.EWAS.by","Av.2step.EWAS_LM.none","Av.2step.EWAS_LM.bon",
                                                        "Av.2step.EWAS_LM.bh","Av.2step.EWAS_LM.by","Av.2step.ENET_MIN","Av.2step.ENET_OPT","Av.2step.SPLS","Av.2step.MMPC","Av.2step.DSA"),]


ggplot(data = plot1.av.2step)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1 - Performance to identify the true exposure whatever the true time point - Averaged exposure levels, 2nd step")


# DLNM

data.dlnm<-result.data1.dlnm.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.dlnm$Method)

data.dlnm$color.perf<-as.factor(ifelse(data.dlnm$sensitivity_mean>80 & data.dlnm$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                       ifelse(data.dlnm$sensitivity_mean>70 & data.dlnm$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                              ifelse(data.dlnm$sensitivity_mean>60 & data.dlnm$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot1.dlnm<-data.dlnm[data.dlnm$Method%in%c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                            "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                            "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by"),]


ggplot(data = plot1.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1 - Performance to identify the true exposure whatever the true time point - DLNM")







######################################
###             data2              ###
######################################


# data2 = all time points - raw data

load(file="dataY2andX/exp3/result2.perf.expo.exp3.2023.RData")
result2.perf.expo.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/result2.perf.expo.exp5.2023.RData")
result2.perf.expo.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/result2.perf.expo.exp10.2023.RData")
result2.perf.expo.exp10$nb.pred<-"10"

result.data2.100<-rbind(result2.perf.expo.exp3,result2.perf.expo.exp5,result2.perf.expo.exp10)
colnames(result.data2.100)

result.data2.100 <- result.data2.100[result.data2.100$Method %in% c("EWAS.TP.none_sum", "EWAS.TP.bon_sum", "EWAS.TP.bh_sum", "EWAS.TP.by_sum", 
                                                                      "EWAS_LM.TP.none_sum", "EWAS_LM.TP.bon_sum", "EWAS_LM.TP.bh_sum", "EWAS_LM.TP.by_sum",
                                                                      "ENET_MIN.TP_sum","ENET_OPT.TP_sum","SPLS_MIN.TP_sum","sNPLS.TP_sum","MMPC.TP_sum","DSA.TP_sum"),]

result.data2.100$Method <- factor(result.data2.100$Method, levels = c("EWAS.TP.none_sum", "EWAS.TP.bon_sum", "EWAS.TP.bh_sum", "EWAS.TP.by_sum", 
                                                                      "EWAS_LM.TP.none_sum", "EWAS_LM.TP.bon_sum", "EWAS_LM.TP.bh_sum", "EWAS_LM.TP.by_sum",
                                                                      "ENET_MIN.TP_sum","ENET_OPT.TP_sum","SPLS_MIN.TP_sum","sNPLS.TP_sum","MMPC.TP_sum","DSA.TP_sum"))

levels(result.data2.100$Method) <- list(Raw.EWAS.none="EWAS.TP.none_sum", Raw.EWAS.bon="EWAS.TP.bon_sum", Raw.EWAS.bh="EWAS.TP.bh_sum", Raw.EWAS.by="EWAS.TP.by_sum", 
                                        Raw.EWAS_LM.none="EWAS_LM.TP.none_sum", Raw.EWAS_LM.bon="EWAS_LM.TP.bon_sum", Raw.EWAS_LM.bh="EWAS_LM.TP.bh_sum", Raw.EWAS_LM.by="EWAS_LM.TP.by_sum",
                                        Raw.ENET_MIN="ENET_MIN.TP_sum",Raw.ENET_OPT="ENET_OPT.TP_sum",Raw.SPLS="SPLS_MIN.TP_sum",Raw.sNPLS="sNPLS.TP_sum",Raw.MMPC="MMPC.TP_sum",Raw.DSA="DSA.TP_sum")

result.data2.100$nb.pred <- factor(result.data2.100$nb.pred, levels = c("3","5","10"))


# data2 = all time points - average (=mean variable)

#### Mistake (the following datasets corresponds to the performance to identify the true expo at the 1-step of the 2-step approach applied with all methods)
load(file="dataY2andX/exp3/twostep/result.data2av.exp3.100.RData")
result.data2av.exp3.100$nb.pred<-"3"
load(file="dataY2andX/exp5/twostep/result.data2av.exp5.100.RData")
result.data2av.exp5.100$nb.pred<-"5"
load(file="dataY2andX/exp10/twostep/result.data2av.exp10.100.RData")
result.data2av.exp10.100$nb.pred<-"10"

result.data2av.100<-rbind(result.data2av.exp3.100,result.data2av.exp5.100,result.data2av.exp10.100)


colnames(result.data2av.100)
unique(result.data2av.100$Method)

result.data2av.100$Method <- factor(result.data2av.100$Method, levels = c("EWAS.TP.none", "EWAS.TP.bon", "EWAS.TP.bh", 
                                                                          "EWAS.TP.by","EWAS_LM.TP.none", "EWAS_LM.TP.bon",
                                                                          "EWAS_LM.TP.bh", "EWAS_LM.TP.by", "ENET_MIN.TP",
                                                                          "ENET_OPT.TP","SPLS_MIN.TP","MMPC.TP","DSA.TP"))

levels(result.data2av.100$Method) <- list(Av.EWAS.none="EWAS.TP.none", Av.EWAS.bon="EWAS.TP.bon", Av.EWAS.bh="EWAS.TP.bh", Av.EWAS.by="EWAS.TP.by", 
                                          Av.EWAS_LM.none="EWAS_LM.TP.none", Av.EWAS_LM.bon="EWAS_LM.TP.bon", Av.EWAS_LM.bh="EWAS_LM.TP.bh", Av.EWAS_LM.by="EWAS_LM.TP.by",
                                          Av.ENET_MIN="ENET_MIN.TP",Av.ENET_OPT="ENET_OPT.TP",Av.SPLS="SPLS_MIN.TP",Av.MMPC="MMPC.TP",Av.DSA="DSA.TP")


result.data2av.100$nb.pred <- factor(result.data2av.100$nb.pred, levels = c("3","5","10"))



# data2 = all time points - average - results obtained after the 2nd step (=raw_variable)

#### the following datasets corresponds to the performance to identify the true expo after the 2-step approach
load(file="dataY2andX/exp3/twostep/result2.perf.expo.av.exp3.RData")
result2.perf.expo.av.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/twostep/result2.perf.expo.av.exp5.RData")
result2.perf.expo.av.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/twostep/result2.perf.expo.av.exp10.RData")
result2.perf.expo.av.exp10$nb.pred<-"10"

result.data2av.2step.100<-rbind(result2.perf.expo.av.exp3,result2.perf.expo.av.exp5,result2.perf.expo.av.exp10)


colnames(result.data2av.2step.100)
unique(result.data2av.2step.100$Method)

result.data2av.2step.100$Method <- factor(result.data2av.2step.100$Method, levels = c("EWAS.TP.none_sum", "EWAS.TP.bon_sum", "EWAS.TP.bh_sum", 
                                                                                      "EWAS.TP.by_sum","EWAS_LM.TP.none_sum", "EWAS_LM.TP.bon_sum",
                                                                                      "EWAS_LM.TP.bh_sum", "EWAS_LM.TP.by_sum", "ENET_MIN.TP_sum",
                                                                                      "ENET_OPT.TP_sum","SPLS_MIN.TP_sum","MMPC.TP_sum","DSA.TP_sum"))

levels(result.data2av.2step.100$Method) <- list(Av.2step.EWAS.none="EWAS.TP.none_sum", Av.2step.EWAS.bon="EWAS.TP.bon_sum", Av.2step.EWAS.bh="EWAS.TP.bh_sum", Av.2step.EWAS.by="EWAS.TP.by_sum", 
                                                Av.2step.EWAS_LM.none="EWAS_LM.TP.none_sum", Av.2step.EWAS_LM.bon="EWAS_LM.TP.bon_sum", Av.2step.EWAS_LM.bh="EWAS_LM.TP.bh_sum", Av.2step.EWAS_LM.by="EWAS_LM.TP.by_sum",
                                                Av.2step.ENET_MIN="ENET_MIN.TP_sum",Av.2step.ENET_OPT="ENET_OPT.TP_sum",Av.2step.SPLS="SPLS_MIN.TP_sum",Av.2step.MMPC="MMPC.TP_sum",Av.2step.DSA="DSA.TP_sum")


result.data2av.2step.100$nb.pred <- factor(result.data2av.2step.100$nb.pred, levels = c("3","5","10"))




# data2 - all time point - DLNM - 1st step (crossbasis)

load(file="dataY2andX/exp3/dlnm/result2.perf.expo.dlnm.exp3.100.RData")
result2.perf.expo.dlnm.exp3.100$nb.pred<-"3"
load(file="dataY2andX/exp5/dlnm/result2.perf.expo.dlnm.exp5.100.RData")
result2.perf.expo.dlnm.exp5.100$nb.pred<-"5"
load(file="dataY2andX/exp10/dlnm/result2.perf.expo.dlnm.exp10.100.RData")
result2.perf.expo.dlnm.exp10.100$nb.pred<-"10"

result.data2.dlnm.100<-rbind(result2.perf.expo.dlnm.exp3.100,result2.perf.expo.dlnm.exp5.100,result2.perf.expo.dlnm.exp10.100)

colnames(result.data2.dlnm.100)
unique(result.data2.dlnm.100$Method)

result.data2.dlnm.100$Method <- factor(result.data2.dlnm.100$Method, levels = c("DLNMpen.none.TP01.100_sum","DLNMpen.bonf.TP01.100_sum","DLNMpen.bh.TP01.100_sum","DLNMpen.by.TP01.100_sum",
                                                                                "DLNMselect.none.TP01.100_sum","DLNMselect.bonf.TP01.100_sum","DLNMselect.bh.TP01.100_sum","DLNMselect.by.TP01.100_sum",
                                                                                "DLNMselectback.none.TP01.100_sum","DLNMselectback.bonf.TP01.100_sum","DLNMselectback.bh.TP01.100_sum","DLNMselectback.by.TP01.100_sum",
                                                                                "DLNMselectforward.none.TP01.100_sum","DLNMselectforward.bonf.TP01.100_sum","DLNMselectforward.bh.TP01.100_sum","DLNMselectforward.by.TP01.100_sum"))

levels(result.data2.dlnm.100$Method) <- list(DLNMpen.none="DLNMpen.none.TP01.100_sum", DLNMpen.bonf="DLNMpen.bonf.TP01.100_sum", DLNMpen.bh="DLNMpen.bh.TP01.100_sum", DLNMpen.by="DLNMpen.by.TP01.100_sum", 
                                             DLNMselect.none="DLNMselect.none.TP01.100_sum", DLNMselect.bonf="DLNMselect.bonf.TP01.100_sum", DLNMselect.bh="DLNMselect.bh.TP01.100_sum", DLNMselect.by="DLNMselect.by.TP01.100_sum",
                                             DLNMselectback.none="DLNMselectback.none.TP01.100_sum",DLNMselectback.bonf="DLNMselectback.bonf.TP01.100_sum",DLNMselectback.bh="DLNMselectback.bh.TP01.100_sum",DLNMselectback.by="DLNMselectback.by.TP01.100_sum",
                                             DLNMselectforward.none="DLNMselectforward.none.TP01.100_sum",DLNMselectforward.bonf="DLNMselectforward.bonf.TP01.100_sum",DLNMselectforward.bh="DLNMselectforward.bh.TP01.100_sum",DLNMselectforward.by="DLNMselectforward.by.TP01.100_sum")

result.data2.dlnm.100$nb.pred <- factor(result.data2.dlnm.100$nb.pred, levels = c("3","5","10"))

result.data2.dlnm.100<-result.data2.dlnm.100[! result.data2.dlnm.100$Method %in% c("DLNMselectforward.bonf", "DLNMselectforward.bh", "DLNMselectforward.by"),]

# data2 - all time point - DLNM - 2nd step (lag)

load(file="dataY2andX/exp3/dlnm/result2.perf.expo.dlnm.exp3.RData")
result2.perf.expo.dlnm.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/dlnm/result2.perf.expo.dlnm.exp5.RData")
result2.perf.expo.dlnm.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/dlnm/result2.perf.expo.dlnm.exp10.RData")
result2.perf.expo.dlnm.exp10$nb.pred<-"10"

result.data2.dlnm.step2.100<-rbind(result2.perf.expo.dlnm.exp3,result2.perf.expo.dlnm.exp5,result2.perf.expo.dlnm.exp10)

colnames(result.data2.dlnm.step2.100)
unique(result.data2.dlnm.step2.100$Method)

result.data2.dlnm.step2.100$Method <- factor(result.data2.dlnm.step2.100$Method, levels = c("DLNMpen.none_sum","DLNMpen.bonf_sum","DLNMpen.bh_sum","DLNMpen.by_sum",
                                                                                            "DLNMselect.none_sum","DLNMselect.bonf_sum","DLNMselect.bh_sum","DLNMselect.by_sum",
                                                                                            "DLNMselectback.none_sum","DLNMselectback.bonf_sum","DLNMselectback.bh_sum","DLNMselectback.by_sum",
                                                                                            "DLNMselectforward.none_sum","DLNMselectforward.bonf_sum","DLNMselectforward.bh_sum","DLNMselectforward.by_sum"))

levels(result.data2.dlnm.step2.100$Method) <- list(DLNMpen.2step.none="DLNMpen.none_sum", DLNMpen.2step.bonf="DLNMpen.bonf_sum", DLNMpen.2step.bh="DLNMpen.bh_sum", DLNMpen.2step.by="DLNMpen.by_sum", 
                                                   DLNMselect.2step.none="DLNMselect.none_sum", DLNMselect.2step.bonf="DLNMselect.bonf_sum", DLNMselect.2step.bh="DLNMselect.bh_sum", DLNMselect.2step.by="DLNMselect.by_sum",
                                                   DLNMselectback.2step.none="DLNMselectback.none_sum",DLNMselectback.2step.bonf="DLNMselectback.bonf_sum",DLNMselectback.2step.bh="DLNMselectback.bh_sum",DLNMselectback.2step.by="DLNMselectback.by_sum",
                                                   DLNMselectforward.2step.none="DLNMselectforward.none_sum",DLNMselectforward.2step.bonf="DLNMselectforward.bonf_sum",DLNMselectforward.2step.bh="DLNMselectforward.bh_sum",DLNMselectforward.2step.by="DLNMselectforward.by_sum")

result.data2.dlnm.step2.100$nb.pred <- factor(result.data2.dlnm.step2.100$nb.pred, levels = c("3","5","10"))

result.data2.dlnm.step2.100<-result.data2.dlnm.step2.100[! result.data2.dlnm.step2.100$Method %in% c("DLNMselectforward.2step.bonf", "DLNMselectforward.2step.bh", "DLNMselectforward.2step.by"),]


# data2 - all time point - DLNM AVG - 1st step (crossbasis)

load(file="dataY2andX/exp3/dlnm/result2.perf.expo.dlnm.AVG.exp3.100.RData")
result2.perf.expo.dlnm.AVG.exp3.100$nb.pred<-"3"
load(file="dataY2andX/exp5/dlnm/result2.perf.expo.dlnm.AVG.exp5.100.RData")
result2.perf.expo.dlnm.AVG.exp5.100$nb.pred<-"5"
load(file="dataY2andX/exp10/dlnm/result2.perf.expo.dlnm.AVG.exp10.100.RData")
result2.perf.expo.dlnm.AVG.exp10.100$nb.pred<-"10"

result.data2.dlnm.AVG.100<-rbind(result2.perf.expo.dlnm.AVG.exp3.100,result2.perf.expo.dlnm.AVG.exp5.100,result2.perf.expo.dlnm.AVG.exp10.100)

colnames(result.data2.dlnm.AVG.100)
unique(result.data2.dlnm.AVG.100$Method)

result.data2.dlnm.AVG.100$Method <- factor(result.data2.dlnm.AVG.100$Method, levels = c("DLNMpen.TP01.100_sum","DLNMselect.TP01.100_sum","DLNMselectback.TP01.100_sum"))

levels(result.data2.dlnm.AVG.100$Method) <- list(DLNMpen.AVG="DLNMpen.TP01.100_sum", DLNMselect.AVG="DLNMselect.TP01.100_sum", DLNMselectback.AVG="DLNMselectback.TP01.100_sum")

result.data2.dlnm.AVG.100$nb.pred <- factor(result.data2.dlnm.AVG.100$nb.pred, levels = c("3","5","10"))


# data2 - all time point - DLNM AVG - 2nd step (lag)

load(file="dataY2andX/exp3/dlnm/result2.perf.expo.dlnm.AVG.exp3.RData")
result2.perf.expo.dlnm.AVG.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/dlnm/result2.perf.expo.dlnm.AVG.exp5.RData")
result2.perf.expo.dlnm.AVG.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/dlnm/result2.perf.expo.dlnm.AVG.exp10.RData")
result2.perf.expo.dlnm.AVG.exp10$nb.pred<-"10"

result.data2.dlnm.AVG.step2.100<-rbind(result2.perf.expo.dlnm.AVG.exp3,result2.perf.expo.dlnm.AVG.exp5,result2.perf.expo.dlnm.AVG.exp10)

colnames(result.data2.dlnm.AVG.step2.100)
unique(result.data2.dlnm.AVG.step2.100$Method)

result.data2.dlnm.AVG.step2.100$Method <- factor(result.data2.dlnm.AVG.step2.100$Method, levels = c("DLNMpen.none_sum","DLNMpen.bonf_sum","DLNMpen.bh_sum","DLNMpen.by_sum",
                                                                                                    "DLNMselect.none_sum","DLNMselect.bonf_sum","DLNMselect.bh_sum","DLNMselect.by_sum",
                                                                                                    "DLNMselectback.none_sum","DLNMselectback.bonf_sum","DLNMselectback.bh_sum","DLNMselectback.by_sum"))

levels(result.data2.dlnm.AVG.step2.100$Method) <- list(DLNMpen.AVG.2step.none="DLNMpen.none_sum", DLNMpen.AVG.2step.bonf="DLNMpen.bonf_sum", DLNMpen.AVG.2step.bh="DLNMpen.bh_sum", DLNMpen.AVG.2step.by="DLNMpen.by_sum", 
                                                       DLNMselect.AVG.2step.none="DLNMselect.none_sum", DLNMselect.AVG.2step.bonf="DLNMselect.bonf_sum", DLNMselect.AVG.2step.bh="DLNMselect.bh_sum", DLNMselect.AVG.2step.by="DLNMselect.by_sum",
                                                       DLNMselectback.AVG.2step.none="DLNMselectback.none_sum",DLNMselectback.AVG.2step.bonf="DLNMselectback.bonf_sum",DLNMselectback.AVG.2step.bh="DLNMselectback.bh_sum",DLNMselectback.AVG.2step.by="DLNMselectback.by_sum")

result.data2.dlnm.AVG.step2.100$nb.pred <- factor(result.data2.dlnm.AVG.step2.100$nb.pred, levels = c("3","5","10"))



# total
data2.all.100<-rbind(result.data2.100,result.data2av.100,result.data2av.2step.100,result.data2.dlnm.100,result.data2.dlnm.step2.100,result.data2.dlnm.AVG.100,result.data2.dlnm.AVG.step2.100)


data2.100<-data2.all.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity),
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate)
  )


data2.100$Sensitivity<-paste0(round(data2.100$sensitivity_mean,1)," (",round(data2.100$sensitivity_sd,1),")")
data2.100$FDR<-paste0(round(data2.100$false.disc.rate_mean,1)," (",round(data2.100$false.disc.rate_sd,1),")")

save(data2.100,file="Summary/data2.100.2023.Rdata")

library(reshape2)
summary_data2.100<-dcast(melt(data2.100, id.vars=c("Method", "nb.pred")), Method~variable+nb.pred)

write.csv2(summary_data2.100,file="Summary/summary.data2.100.2023.csv")
save(summary_data2.100,file="Summary/summary.data2.100.2023.Rdata")



################# PLOT


data<-result.data2.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data$Method)


data$color.perf<-as.factor(ifelse(data$sensitivity_mean>80 & data$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                  ifelse(data$sensitivity_mean>70 & data$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                         ifelse(data$sensitivity_mean>60 & data$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))

# raw data
plot2.raw<-data[data$Method%in%c("Raw.EWAS.none","Raw.EWAS.bon","Raw.EWAS.bh","Raw.EWAS.by","Raw.EWAS_LM.none","Raw.EWAS_LM.bon",
                                 "Raw.EWAS_LM.bh","Raw.EWAS_LM.by","Raw.ENET_MIN","Raw.ENET_OPT","Raw.SPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

ggplot(data = plot2.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2: A single time point associated with Y - Raw data - Performance to identify the true exposure whatever the true time point")


# averaged levels - 1st step

data.av<-result.data2av.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.av$Method)

data.av$color.perf<-as.factor(ifelse(data.av$sensitivity_mean>80 & data.av$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                     ifelse(data.av$sensitivity_mean>70 & data.av$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                            ifelse(data.av$sensitivity_mean>60 & data.av$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot2.av<-data.av[data.av$Method%in%c("Av.EWAS.none","Av.EWAS.bon","Av.EWAS.bh","Av.EWAS.by","Av.EWAS_LM.none","Av.EWAS_LM.bon",
                                      "Av.EWAS_LM.bh","Av.EWAS_LM.by","Av.ENET_MIN","Av.ENET_OPT","Av.SPLS","Av.MMPC","Av.DSA"),]


ggplot(data = plot2.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2 - Performance to identify the true exposure whatever the true time point - Averaged exposure levels, 1st step")



# averaged levels - 2-step

data.av.2step<-result.data2av.2step.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.av.2step$Method)

data.av.2step$color.perf<-as.factor(ifelse(data.av.2step$sensitivity_mean>80 & data.av.2step$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                           ifelse(data.av.2step$sensitivity_mean>70 & data.av.2step$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                                  ifelse(data.av.2step$sensitivity_mean>60 & data.av.2step$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot2.av.2step<-data.av.2step[data.av.2step$Method%in%c("Av.2step.EWAS.none","Av.2step.EWAS.bon","Av.2step.EWAS.bh","Av.2step.EWAS.by","Av.2step.EWAS_LM.none","Av.2step.EWAS_LM.bon",
                                                        "Av.2step.EWAS_LM.bh","Av.2step.EWAS_LM.by","Av.2step.ENET_MIN","Av.2step.ENET_OPT","Av.2step.SPLS","Av.2step.MMPC","Av.2step.DSA"),]


ggplot(data = plot2.av.2step)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2 - Performance to identify the true exposure whatever the true time point - Averaged exposure levels, 2nd step")



# DLNM

data.dlnm<-result.data2.dlnm.100 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.dlnm$Method)

data.dlnm$color.perf<-as.factor(ifelse(data.dlnm$sensitivity_mean>80 & data.dlnm$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                       ifelse(data.dlnm$sensitivity_mean>70 & data.dlnm$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                              ifelse(data.dlnm$sensitivity_mean>60 & data.dlnm$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot2.dlnm<-data.dlnm[data.dlnm$Method%in%c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                            "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                            "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by"),]


ggplot(data = plot2.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2 - Performance to identify the true exposure whatever the true time point - DLNM")












################################### Performance to identify the true exposure at the true time point ####################################





######################################
###             DATA1              ###
######################################

# Data1 = all time points - raw data

load(file="dataY1andX/exp3/result.data1.exp3.2023.RData")
result.data1.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/result.data1.exp5.2023.RData")
result.data1.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/result.data1.exp10.2023.RData")
result.data1.exp10$nb.pred<-"10"

result.data1<-rbind(result.data1.exp3,result.data1.exp5,result.data1.exp10)
colnames(result.data1)
unique(result.data1$Method)

result.data1$Method <- factor(result.data1$Method, levels = c("EWAS.TP.none", "EWAS.TP.bon", "EWAS.TP.bh", "EWAS.TP.by", 
                                                              "EWAS_LM.TP.none", "EWAS_LM.TP.bon", "EWAS_LM.TP.bh", "EWAS_LM.TP.by",
                                                              "ENET_MIN.TP","ENET_OPT.TP","SPLS_MIN.TP","sNPLS.TP","MMPC.TP","DSA.TP",
                                                              "DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                              "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                              "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                              "DLNMselectforward.none","DLNMselectforward.bonf", "DLNMselectforward.bh", "DLNMselectforward.by"))

levels(result.data1$Method) <- list(Raw.ExWAS.none="EWAS.TP.none", Raw.Exwas.bon="EWAS.TP.bon", Raw.ExWAS.bh="EWAS.TP.bh", Raw.ExWAS.by="EWAS.TP.by", 
                                    Raw.ExWAS.MLM.none="EWAS_LM.TP.none", Raw.Exwas.MLM.bon="EWAS_LM.TP.bon", Raw.ExWAS.MLM.bh="EWAS_LM.TP.bh", Raw.ExWAS.MLM.by="EWAS_LM.TP.by",
                                    Raw.ENET.min="ENET_MIN.TP",Raw.ENET.opt="ENET_OPT.TP",Raw.sPLS="SPLS_MIN.TP",Raw.sNPLS="sNPLS.TP",Raw.MMPC="MMPC.TP",Raw.DSA="DSA.TP",
                                    DLNMpen.none="DLNMpen.none", DLNMpen.bonf="DLNMpen.bonf", DLNMpen.bh="DLNMpen.bh", DLNMpen.by="DLNMpen.by", 
                                    DLNMselect.none="DLNMselect.none", DLNMselect.bonf="DLNMselect.bonf", DLNMselect.bh="DLNMselect.bh", DLNMselect.by="DLNMselect.by",
                                    DLNMselectback.none="DLNMselectback.none",DLNMselectback.bonf="DLNMselectback.bonf",DLNMselectback.bh="DLNMselectback.bh",DLNMselectback.by="DLNMselectback.by",
                                    DLNMselectforward.none="DLNMselectforward.none",DLNMselectforward.bonf="DLNMselectforward.bonf",DLNMselectforward.bh="DLNMselectforward.bh",DLNMselectforward.by="DLNMselectforward.by")

result.data1$nb.pred <- factor(result.data1$nb.pred, levels = c("3","5","10"))



# Data1 = all time points - average reduced

load(file="dataY1andX/exp3/twostep/result.data1.red.exp3.RData")
result.data1.red.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/twostep/result.data1.red.exp5.RData")
result.data1.red.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/twostep/result.data1.red.exp10.RData")
result.data1.red.exp10$nb.pred<-"10"


result.data1av<-rbind(result.data1.red.exp3,result.data1.red.exp5,result.data1.red.exp10)
colnames(result.data1av)
unique(result.data1av$Method)

result.data1av$Method <- factor(result.data1av$Method, levels = c("EWAS.TP.none", "EWAS.TP.bon", "EWAS.TP.bh", "EWAS.TP.by", 
                                                                  "EWAS_LM.TP.none", "EWAS_LM.TP.bon", "EWAS_LM.TP.bh", "EWAS_LM.TP.by",
                                                                  "ENET_MIN.TP","ENET_OPT.TP","SPLS_MIN.TP","MMPC.TP","DSA.TP"))

levels(result.data1av$Method) <- list(Av.ExWAS.none="EWAS.TP.none", Av.Exwas.bon="EWAS.TP.bon", Av.ExWAS.bh="EWAS.TP.bh", Av.ExWAS.by="EWAS.TP.by", 
                                      Av.ExWAS.MLM.none="EWAS_LM.TP.none", Av.Exwas.MLM.bon="EWAS_LM.TP.bon", Av.ExWAS.MLM.bh="EWAS_LM.TP.bh", Av.ExWAS.MLM.by="EWAS_LM.TP.by",
                                      Av.ENET.min="ENET_MIN.TP",Av.ENET.opt="ENET_OPT.TP",Av.sPLS="SPLS_MIN.TP",Av.MMPC="MMPC.TP",Av.DSA="DSA.TP")

result.data1av$nb.pred <- factor(result.data1av$nb.pred, levels = c("3","5","10"))

result.data1av<-result.data1av[,-which(colnames(result.data1av)%in%c("n.candidate.red"))]



####  Data1 = all time points - DLNM

# load(file="dataY1andX/exp3/dlnm/result.data1.dlnm.exp3.RData")
# result.data1.dlnm.exp3$nb.pred<-"3"
# load(file="dataY1andX/exp5/dlnm/result.data1.dlnm.exp5.RData")
# result.data1.dlnm.exp5$nb.pred<-"5"
# load(file="dataY1andX/exp10/dlnm/result.data1.dlnm.exp10.RData")
# result.data1.dlnm.exp10$nb.pred<-"10"
# 
# 
# result.data1.dlnm<-rbind(result.data1.dlnm.exp3,result.data1.dlnm.exp5,result.data1.dlnm.exp10)
# colnames(result.data1.dlnm)
# unique(result.data1.dlnm$Method)
# 
# result.data1.dlnm$Method <- factor(result.data1.dlnm$Method, levels = c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
#                                                                         "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
#                                                                         "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
#                                                                         "DLNMselectforward.none","DLNMselectforward.bonf", "DLNMselectforward.bh", "DLNMselectforward.by"))
# 
# levels(result.data1.dlnm$Method) <- list(DLNMpen.none="DLNMpen.none", DLNMpen.bonf="DLNMpen.bonf", DLNMpen.bh="DLNMpen.bh", DLNMpen.by="DLNMpen.by", 
#                                          DLNMselect.none="DLNMselect.none", DLNMselect.bonf="DLNMselect.bonf", DLNMselect.bh="DLNMselect.bh", DLNMselect.by="DLNMselect.by",
#                                          DLNMselectback.none="DLNMselectback.none",DLNMselectback.bonf="DLNMselectback.bonf",DLNMselectback.bh="DLNMselectback.bh",DLNMselectback.by="DLNMselectback.by",
#                                          DLNMselectforward.none="DLNMselectforward.none",DLNMselectforward.bonf="DLNMselectforward.bonf",DLNMselectforward.bh="DLNMselectforward.bh",DLNMselectforward.by="DLNMselectforward.by")
# 
# result.data1.dlnm$nb.pred <- factor(result.data1.dlnm$nb.pred, levels = c("3","5","10"))
# 



####  Data1 = all time points - DLNM AVG

load(file="dataY1andX/exp3/dlnm/result.data1.dlnm.AVG.exp3.RData")
result.data1.dlnm.AVG.exp3$nb.pred<-"3"
load(file="dataY1andX/exp5/dlnm/result.data1.dlnm.AVG.exp5.RData")
result.data1.dlnm.AVG.exp5$nb.pred<-"5"
load(file="dataY1andX/exp10/dlnm/result.data1.dlnm.AVG.exp10.RData")
result.data1.dlnm.AVG.exp10$nb.pred<-"10"


result.data1.dlnm.AVG<-rbind(result.data1.dlnm.AVG.exp3,result.data1.dlnm.AVG.exp5,result.data1.dlnm.AVG.exp10)
colnames(result.data1.dlnm.AVG)
unique(result.data1.dlnm.AVG$Method)

result.data1.dlnm.AVG$Method <- factor(result.data1.dlnm.AVG$Method, levels = c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                                                "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                                                "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by"))

levels(result.data1.dlnm.AVG$Method) <- list(DLNMpen.AVG.none="DLNMpen.none", DLNMpen.AVG.bonf="DLNMpen.bonf", DLNMpen.AVG.bh="DLNMpen.bh", DLNMpen.AVG.by="DLNMpen.by", 
                                             DLNMselect.AVG.none="DLNMselect.none", DLNMselect.AVG.bonf="DLNMselect.bonf", DLNMselect.AVG.bh="DLNMselect.bh", DLNMselect.AVG.by="DLNMselect.by",
                                             DLNMselectback.AVG.none="DLNMselectback.none",DLNMselectback.AVG.bonf="DLNMselectback.bonf",DLNMselectback.AVG.bh="DLNMselectback.bh",DLNMselectback.AVG.by="DLNMselectback.by")

result.data1.dlnm.AVG$nb.pred <- factor(result.data1.dlnm.AVG$nb.pred, levels = c("3","5","10"))



### total


data1.all<-rbind(result.data1,result.data1av,result.data1.dlnm.AVG)

save(data1.all,file="Summary/data1.all.2023.Rdata")

data<-data1.all %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity),
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate)  
    )


data$Sensitivity<-paste0(round(data$sensitivity_mean,1)," (",round(data$sensitivity_sd,1),")")
data$FDR<-paste0(round(data$false.disc.rate_mean,1)," (",round(data$false.disc.rate_sd,1),")")


library(reshape2)
summary_data1<-dcast(melt(data, id.vars=c("Method", "nb.pred")), Method~variable+nb.pred)


write.csv2(summary_data1,file="Summary/summary.data1.2023.csv")
save(summary_data1,file="Summary/summary.data1.2023.Rdata")



############ PLOT

data<-result.data1 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data$Method)


data$color.perf<-as.factor(ifelse(data$sensitivity_mean>80 & data$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                  ifelse(data$sensitivity_mean>70 & data$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                         ifelse(data$sensitivity_mean>60 & data$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))

# raw data
plot1.raw<-data[data$Method%in%c("Raw.ExWAS.none","Raw.Exwas.bon","Raw.ExWAS.bh","Raw.ExWAS.by","Raw.ExWAS.MLM.none","Raw.Exwas.MLM.bon",
                                 "Raw.ExWAS.MLM.bh","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

ggplot(data = plot1.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1: All time points associated with Y - Raw data - Performance to identify the true exposure at the true time point")



# averaged levels

data.av<-result.data1av %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.av$Method)

data.av$color.perf<-as.factor(ifelse(data.av$sensitivity_mean>80 & data.av$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                     ifelse(data.av$sensitivity_mean>70 & data.av$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                            ifelse(data.av$sensitivity_mean>60 & data.av$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot1.av<-data.av[data.av$Method%in%c("Av.ExWAS.none","Av.Exwas.bon","Av.ExWAS.bh","Av.ExWAS.by","Av.ExWAS.MLM.none","Av.Exwas.MLM.bon",
                                      "Av.ExWAS.MLM.bh","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA"),]


ggplot(data = plot1.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1 - Performance to identify the true exposure at the true time point - Averaged exposure levels")



# DLNM

data.dlnm<-result.data1 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.dlnm$Method)

data.dlnm$color.perf<-as.factor(ifelse(data.dlnm$sensitivity_mean>80 & data.dlnm$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                       ifelse(data.dlnm$sensitivity_mean>70 & data.dlnm$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                              ifelse(data.dlnm$sensitivity_mean>60 & data.dlnm$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot1.dlnm<-data.dlnm[data.dlnm$Method%in%c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                            "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                            "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by"),]


ggplot(data = plot1.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 1 - Performance to identify the true exposure at the true time point - DLNM")





######################################
###             data2              ###
######################################

# data2 = all time points - raw data

load(file="dataY2andX/exp3/result.data2.exp3.2023.RData")
result.data2.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/result.data2.exp5.2023.RData")
result.data2.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/result.data2.exp10.2023.RData")
result.data2.exp10$nb.pred<-"10"

result.data2<-rbind(result.data2.exp3,result.data2.exp5,result.data2.exp10)
colnames(result.data2)
unique(result.data2$Method)

result.data2$Method <- factor(result.data2$Method, levels = c("EWAS.TP.none", "EWAS.TP.bon", "EWAS.TP.bh", "EWAS.TP.by", 
                                                              "EWAS_LM.TP.none", "EWAS_LM.TP.bon", "EWAS_LM.TP.bh", "EWAS_LM.TP.by",
                                                              "ENET_MIN.TP","ENET_OPT.TP","SPLS_MIN.TP","sNPLS.TP","MMPC.TP","DSA.TP",
                                                              "DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                              "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                              "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                              "DLNMselectforward.none","DLNMselectforward.bonf", "DLNMselectforward.bh", "DLNMselectforward.by"))

levels(result.data2$Method) <- list(Raw.ExWAS.none="EWAS.TP.none", Raw.Exwas.bon="EWAS.TP.bon", Raw.ExWAS.bh="EWAS.TP.bh", Raw.ExWAS.by="EWAS.TP.by", 
                                    Raw.ExWAS.MLM.none="EWAS_LM.TP.none", Raw.Exwas.MLM.bon="EWAS_LM.TP.bon", Raw.ExWAS.MLM.bh="EWAS_LM.TP.bh", Raw.ExWAS.MLM.by="EWAS_LM.TP.by",
                                    Raw.ENET.min="ENET_MIN.TP",Raw.ENET.opt="ENET_OPT.TP",Raw.sPLS="SPLS_MIN.TP",Raw.sNPLS="sNPLS.TP",Raw.MMPC="MMPC.TP",Raw.DSA="DSA.TP",
                                    DLNMpen.none="DLNMpen.none", DLNMpen.bonf="DLNMpen.bonf", DLNMpen.bh="DLNMpen.bh", DLNMpen.by="DLNMpen.by", 
                                    DLNMselect.none="DLNMselect.none", DLNMselect.bonf="DLNMselect.bonf", DLNMselect.bh="DLNMselect.bh", DLNMselect.by="DLNMselect.by",
                                    DLNMselectback.none="DLNMselectback.none",DLNMselectback.bonf="DLNMselectback.bonf",DLNMselectback.bh="DLNMselectback.bh",DLNMselectback.by="DLNMselectback.by",
                                    DLNMselectforward.none="DLNMselectforward.none",DLNMselectforward.bonf="DLNMselectforward.bonf",DLNMselectforward.bh="DLNMselectforward.bh",DLNMselectforward.by="DLNMselectforward.by")

result.data2$nb.pred <- factor(result.data2$nb.pred, levels = c("3","5","10"))



# data2 = all time points - average reduced

load(file="dataY2andX/exp3/twostep/result.data2.red.exp3.RData")
result.data2.red.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/twostep/result.data2.red.exp5.RData")
result.data2.red.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/twostep/result.data2.red.exp10.RData")
result.data2.red.exp10$nb.pred<-"10"


result.data2av<-rbind(result.data2.red.exp3,result.data2.red.exp5,result.data2.red.exp10)
colnames(result.data2av)
unique(result.data2av$Method)

result.data2av$Method <- factor(result.data2av$Method, levels = c("EWAS.TP.none", "EWAS.TP.bon", "EWAS.TP.bh", "EWAS.TP.by", 
                                                                  "EWAS_LM.TP.none", "EWAS_LM.TP.bon", "EWAS_LM.TP.bh", "EWAS_LM.TP.by",
                                                                  "ENET_MIN.TP","ENET_OPT.TP","SPLS_MIN.TP","MMPC.TP","DSA.TP"))

levels(result.data2av$Method) <- list(Av.ExWAS.none="EWAS.TP.none", Av.Exwas.bon="EWAS.TP.bon", Av.ExWAS.bh="EWAS.TP.bh", Av.ExWAS.by="EWAS.TP.by", 
                                      Av.ExWAS.MLM.none="EWAS_LM.TP.none", Av.Exwas.MLM.bon="EWAS_LM.TP.bon", Av.ExWAS.MLM.bh="EWAS_LM.TP.bh", Av.ExWAS.MLM.by="EWAS_LM.TP.by",
                                      Av.ENET.min="ENET_MIN.TP",Av.ENET.opt="ENET_OPT.TP",Av.sPLS="SPLS_MIN.TP",Av.MMPC="MMPC.TP",Av.DSA="DSA.TP")

result.data2av$nb.pred <- factor(result.data2av$nb.pred, levels = c("3","5","10"))

result.data2av<-result.data2av[,-which(colnames(result.data2av)%in%c("n.candidate.red"))]


####  data2 = all time points - DLNM
# 
# load(file="dataY2andX/exp3/dlnm/result.data2.dlnm.exp3.RData")
# result.data2.dlnm.exp3$nb.pred<-"3"
# load(file="dataY2andX/exp5/dlnm/result.data2.dlnm.exp5.RData")
# result.data2.dlnm.exp5$nb.pred<-"5"
# load(file="dataY2andX/exp10/dlnm/result.data2.dlnm.exp10.RData")
# result.data2.dlnm.exp10$nb.pred<-"10"
# 
# 
# result.data2.dlnm<-rbind(result.data2.dlnm.exp3,result.data2.dlnm.exp5,result.data2.dlnm.exp10)
# colnames(result.data2.dlnm)
# unique(result.data2.dlnm$Method)
# 
# result.data2.dlnm$Method <- factor(result.data2.dlnm$Method, levels = c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
#                                                                         "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
#                                                                         "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
#                                                                         "DLNMselectforward.none","DLNMselectforward.bonf", "DLNMselectforward.bh", "DLNMselectforward.by"))
# 
# levels(result.data2.dlnm$Method) <- list(DLNMpen.none="DLNMpen.none", DLNMpen.bonf="DLNMpen.bonf", DLNMpen.bh="DLNMpen.bh", DLNMpen.by="DLNMpen.by", 
#                                          DLNMselect.none="DLNMselect.none", DLNMselect.bonf="DLNMselect.bonf", DLNMselect.bh="DLNMselect.bh", DLNMselect.by="DLNMselect.by",
#                                          DLNMselectback.none="DLNMselectback.none",DLNMselectback.bonf="DLNMselectback.bonf",DLNMselectback.bh="DLNMselectback.bh",DLNMselectback.by="DLNMselectback.by",
#                                          DLNMselectforward.none="DLNMselectforward.none",DLNMselectforward.bonf="DLNMselectforward.bonf",DLNMselectforward.bh="DLNMselectforward.bh",DLNMselectforward.by="DLNMselectforward.by")
# 
# result.data2.dlnm$nb.pred <- factor(result.data2.dlnm$nb.pred, levels = c("3","5","10"))
# 



####  data2 = all time points - DLNM AVG

load(file="dataY2andX/exp3/dlnm/result.data2.dlnm.AVG.exp3.RData")
result.data2.dlnm.AVG.exp3$nb.pred<-"3"
load(file="dataY2andX/exp5/dlnm/result.data2.dlnm.AVG.exp5.RData")
result.data2.dlnm.AVG.exp5$nb.pred<-"5"
load(file="dataY2andX/exp10/dlnm/result.data2.dlnm.AVG.exp10.RData")
result.data2.dlnm.AVG.exp10$nb.pred<-"10"


result.data2.dlnm.AVG<-rbind(result.data2.dlnm.AVG.exp3,result.data2.dlnm.AVG.exp5,result.data2.dlnm.AVG.exp10)
colnames(result.data2.dlnm.AVG)
unique(result.data2.dlnm.AVG$Method)

result.data2.dlnm.AVG$Method <- factor(result.data2.dlnm.AVG$Method, levels = c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                                                "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                                                "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by"))

levels(result.data2.dlnm.AVG$Method) <- list(DLNMpen.AVG.none="DLNMpen.none", DLNMpen.AVG.bonf="DLNMpen.bonf", DLNMpen.AVG.bh="DLNMpen.bh", DLNMpen.AVG.by="DLNMpen.by", 
                                             DLNMselect.AVG.none="DLNMselect.none", DLNMselect.AVG.bonf="DLNMselect.bonf", DLNMselect.AVG.bh="DLNMselect.bh", DLNMselect.AVG.by="DLNMselect.by",
                                             DLNMselectback.AVG.none="DLNMselectback.none",DLNMselectback.AVG.bonf="DLNMselectback.bonf",DLNMselectback.AVG.bh="DLNMselectback.bh",DLNMselectback.AVG.by="DLNMselectback.by")

result.data2.dlnm.AVG$nb.pred <- factor(result.data2.dlnm.AVG$nb.pred, levels = c("3","5","10"))



### total


data2.all<-rbind(result.data2,result.data2av,result.data2.dlnm.AVG)

save(data2.all,file="Summary/data2.all.2023.Rdata")

data<-data2.all %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity),
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate)
  )


data$Sensitivity<-paste0(round(data$sensitivity_mean,1)," (",round(data$sensitivity_sd,1),")")
data$FDR<-paste0(round(data$false.disc.rate_mean,1)," (",round(data$false.disc.rate_sd,1),")")


library(reshape2)
summary_data2<-dcast(melt(data, id.vars=c("Method", "nb.pred")), Method~variable+nb.pred)

write.csv2(summary_data2,file="Summary/summary.data2.2023.csv")
save(summary_data2,file="Summary/summary.data2.2023.Rdata")


########### PLOT


data<-result.data2 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data$Method)


data$color.perf<-as.factor(ifelse(data$sensitivity_mean>80 & data$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                  ifelse(data$sensitivity_mean>70 & data$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                         ifelse(data$sensitivity_mean>60 & data$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))

# raw data
plot2.raw<-data[data$Method%in%c("Raw.ExWAS.none","Raw.Exwas.bon","Raw.ExWAS.bh","Raw.ExWAS.by","Raw.ExWAS.MLM.none","Raw.Exwas.MLM.bon",
                                 "Raw.ExWAS.MLM.bh","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

ggplot(data = plot2.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2: A single time point associated with Y - Raw data - Performance to identify the true exposure at the true time point")



# averaged levels

data.av<-result.data2av %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.av$Method)

data.av$color.perf<-as.factor(ifelse(data.av$sensitivity_mean>80 & data.av$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                     ifelse(data.av$sensitivity_mean>70 & data.av$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                            ifelse(data.av$sensitivity_mean>60 & data.av$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot2.av<-data.av[data.av$Method%in%c("Av.ExWAS.none","Av.Exwas.bon","Av.ExWAS.bh","Av.ExWAS.by","Av.ExWAS.MLM.none","Av.Exwas.MLM.bon",
                                      "Av.ExWAS.MLM.bh","Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA"),]


ggplot(data = plot2.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2 - Performance to identify the true exposure at the true time point - Averaged exposure levels")



# DLNM

data.dlnm<-result.data2 %>% 
  group_by(Method,nb.pred) %>% 
  summarise(
    false.disc.rate_mean = mean(false.disc.rate),
    false.disc.rate_sd = sd(false.disc.rate),
    sensitivity_mean = mean(sensitivity),
    sensitivity_sd = sd(sensitivity)
  )

table(data.dlnm$Method)

data.dlnm$color.perf<-as.factor(ifelse(data.dlnm$sensitivity_mean>80 & data.dlnm$false.disc.rate_mean<10,"S>80% and FDR<10%",
                                       ifelse(data.dlnm$sensitivity_mean>70 & data.dlnm$false.disc.rate_mean<20,"S>70% and FDR<20%",
                                              ifelse(data.dlnm$sensitivity_mean>60 & data.dlnm$false.disc.rate_mean<30,"S>60% and FDR<30%","S<60% or FDR>30%"))))


plot2.dlnm<-data.dlnm[data.dlnm$Method%in%c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                            "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                            "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by"),]


ggplot(data = plot2.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf))+
  scale_color_manual(values=c("red","orange","blue","green3"))+
  facet_grid(nb.pred~Method)+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("FDR")+
  ggtitle("Scenario 2 - Performance to identify the true exposure at the true time point - DLNM")


