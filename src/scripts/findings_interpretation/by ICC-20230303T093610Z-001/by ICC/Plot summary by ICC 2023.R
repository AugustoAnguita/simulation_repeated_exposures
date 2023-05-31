library(dplyr)
library(ggplot2)



######### Performance to identify the true exposure at the true time point by ICC

### data1

load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data1_ICC_2023.Rdata")


result.data1.ICC<-result.data1.ICC %>% 
  group_by(Method,expo,ICC) %>% 
  summarise(
    false.disc.rate_mean = mean(na.omit(false.disc.rate)),
    sensitivity_mean = mean(na.omit(sensitivity))
    )


result.data1.ICC$ICC<-as.character(result.data1.ICC$ICC)
result.data1.ICC$ICC<-factor(result.data1.ICC$ICC, levels=c("1","0.1","0.5","0.9"))

result.data1.ICC$color.perf<-as.factor(ifelse(result.data1.ICC$sensitivity_mean>80 & result.data1.ICC$false.disc.rate_mean<20,"S>80% and FDR<20%",
                                              ifelse(result.data1.ICC$sensitivity_mean>70 & result.data1.ICC$false.disc.rate_mean<30,"S>70% and FDR<30%",
                                                     ifelse(result.data1.ICC$sensitivity_mean>60 & result.data1.ICC$false.disc.rate_mean<40,"S>60% and FDR<40%","S<60% or FDR>40%"))))

# raw data

data1.raw<-result.data1.ICC[result.data1.ICC$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh",          
                                                           "Raw.ExWAS.MLM.bh","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt",          
                                                           "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

data1.raw$Method<-factor(data1.raw$Method, levels=c("Raw.ExWAS.none","Raw.Exwas.bon","Raw.ExWAS.bh","Raw.ExWAS.by",          
                                                    "Raw.ExWAS.MLM.none","Raw.Exwas.MLM.bon","Raw.ExWAS.MLM.bh","Raw.ExWAS.MLM.by",
                                                    "Raw.ENET.min","Raw.ENET.opt", "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"))

ggplot(data = data1.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - Raw data \nPerformance to identify the true exposure at the true time point by ICC")



# Av data

data1.av<-result.data1.ICC[result.data1.ICC$Method %in% c( "Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon",          
                                                           "Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh","Av.ExWAS.by",           
                                                           "Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS",               
                                                           "Av.MMPC","Av.DSA"),]

data1.av$Method<-factor(data1.av$Method, levels=c("Av.ExWAS.none","Av.Exwas.bon","Av.ExWAS.bh","Av.ExWAS.by",      
                                                  "Av.ExWAS.MLM.none","Av.Exwas.MLM.bon","Av.ExWAS.MLM.bh","Av.ExWAS.MLM.by",
                                                  "Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA"))


ggplot(data = data1.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - Averaged data \nPerformance to identify the true exposure at the true time point by ICC")


# DLNM

data1.dlnm<-result.data1.ICC[result.data1.ICC$Method %in% c( "DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                             "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                             "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                             "DLNMselectforward.none","DLNMselectforward.bonf"),]

data1.dlnm$Method<-factor(data1.dlnm$Method, levels=c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                      "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                      "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                      "DLNMselectforward.none","DLNMselectforward.bonf"))


ggplot(data = data1.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - DLM \nPerformance to identify the true exposure at the true time point by ICC")


# Av + DLNM

data1.dlnm<-result.data1.ICC[result.data1.ICC$Method %in% c( "Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                             "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                             "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by" ),]

data1.dlnm$Method<-factor(data1.dlnm$Method, levels=c("Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                      "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                      "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by"))


ggplot(data = data1.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - DLM after Av.ExWAS \nPerformance to identify the true exposure at the true time point by ICC")





### data2

load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data2_ICC.2023.Rdata")

result.data2.ICC<-result.data2.ICC %>% 
  group_by(Method,expo,ICC) %>% 
  summarise(
    false.disc.rate_mean = mean(na.omit(false.disc.rate)),
    sensitivity_mean = mean(na.omit(sensitivity))
  )

result.data2.ICC$ICC<-as.character(result.data2.ICC$ICC)
result.data2.ICC$ICC<-factor(result.data2.ICC$ICC, levels=c("1","0.1","0.5","0.9"))

result.data2.ICC$color.perf<-as.factor(ifelse(result.data2.ICC$sensitivity_mean>80 & result.data2.ICC$false.disc.rate_mean<20,"S>80% and FDR<20%",
                                              ifelse(result.data2.ICC$sensitivity_mean>70 & result.data2.ICC$false.disc.rate_mean<30,"S>70% and FDR<30%",
                                                     ifelse(result.data2.ICC$sensitivity_mean>60 & result.data2.ICC$false.disc.rate_mean<40,"S>60% and FDR<40%","S<60% or FDR>40%"))))

# raw data

data2.raw<-result.data2.ICC[result.data2.ICC$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh",          
                                                           "Raw.ExWAS.MLM.bh","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt",          
                                                           "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

data2.raw$Method<-factor(data2.raw$Method, levels=c("Raw.ExWAS.none","Raw.Exwas.bon","Raw.ExWAS.bh","Raw.ExWAS.by",          
                                                    "Raw.ExWAS.MLM.none","Raw.Exwas.MLM.bon","Raw.ExWAS.MLM.bh","Raw.ExWAS.MLM.by",
                                                    "Raw.ENET.min","Raw.ENET.opt", "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"))


ggplot(data = data2.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - Raw data \nPerformance to identify the true exposure at the true time point by ICC")



# Av data

data2.av<-result.data2.ICC[result.data2.ICC$Method %in% c( "Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon",          
                                                           "Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh","Av.ExWAS.by",           
                                                           "Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS",               
                                                           "Av.MMPC","Av.DSA"),]

data2.av$Method<-factor(data2.av$Method, levels=c("Av.ExWAS.none","Av.Exwas.bon","Av.ExWAS.bh","Av.ExWAS.by",      
                                                  "Av.ExWAS.MLM.none","Av.Exwas.MLM.bon","Av.ExWAS.MLM.bh","Av.ExWAS.MLM.by",
                                                  "Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA"))


ggplot(data = data2.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - Averaged data \nPerformance to identify the true exposure at the true time point by ICC")


# DLNM

data2.dlnm<-result.data2.ICC[result.data2.ICC$Method %in% c( "DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                             "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                             "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                             "DLNMselectforward.none","DLNMselectforward.bonf"),]

data2.dlnm$Method<-factor(data2.dlnm$Method, levels=c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                      "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                      "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                      "DLNMselectforward.none","DLNMselectforward.bonf"))


ggplot(data = data2.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - DLM \nPerformance to identify the true exposure at the true time point by ICC")


# Av + DLNM

data2.dlnm<-result.data2.ICC[result.data2.ICC$Method %in% c( "Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                             "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                             "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by" ),]

data2.dlnm$Method<-factor(data2.dlnm$Method, levels=c("Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                      "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                      "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by"))


ggplot(data = data2.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - DLM after Av.ExWAS \nPerformance to identify the true exposure at the true time point by ICC")






######### Performance to identify the true exposure at the true time point by ICC


### data1

load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data1.100_ICC.2023.Rdata")


result.data1.100.ICC<-result.data1.100.ICC %>% 
  group_by(Method,expo,ICC) %>% 
  summarise(
    false.disc.rate_mean = mean(na.omit(false.disc.rate)),
    sensitivity_mean = mean(na.omit(sensitivity))
  )

result.data1.100.ICC$ICC<-as.character(result.data1.100.ICC$ICC)
result.data1.100.ICC$ICC<-factor(result.data1.100.ICC$ICC, levels=c("1","0.1","0.5","0.9"))

result.data1.100.ICC$color.perf<-as.factor(ifelse(result.data1.100.ICC$sensitivity_mean>80 & result.data1.100.ICC$false.disc.rate_mean<20,"S>80% and FDR<20%",
                                                  ifelse(result.data1.100.ICC$sensitivity_mean>70 & result.data1.100.ICC$false.disc.rate_mean<30,"S>70% and FDR<30%",
                                                         ifelse(result.data1.100.ICC$sensitivity_mean>60 & result.data1.100.ICC$false.disc.rate_mean<40,"S>60% and FDR<40%","S<60% or FDR>40%"))))

# raw data

data1.raw<-result.data1.100.ICC[result.data1.100.ICC$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh",          
                                                                   "Raw.ExWAS.MLM.bh","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt",          
                                                                   "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

data1.raw$Method<-factor(data1.raw$Method, levels=c("Raw.ExWAS.none","Raw.Exwas.bon","Raw.ExWAS.bh","Raw.ExWAS.by",          
                                                    "Raw.ExWAS.MLM.none","Raw.Exwas.MLM.bon","Raw.ExWAS.MLM.bh","Raw.ExWAS.MLM.by",
                                                    "Raw.ENET.min","Raw.ENET.opt", "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"))


ggplot(data = data1.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - Raw data \nPerformance to identify the true exposure whatever the time point by ICC")


# Av data

data1.av<-result.data1.100.ICC[result.data1.100.ICC$Method %in% c( "Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon",          
                                                                   "Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh","Av.ExWAS.by",           
                                                                   "Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS",               
                                                                   "Av.MMPC","Av.DSA"),]

data1.av$Method<-factor(data1.av$Method, levels=c("Av.ExWAS.none","Av.Exwas.bon","Av.ExWAS.bh","Av.ExWAS.by",      
                                                  "Av.ExWAS.MLM.none","Av.Exwas.MLM.bon","Av.ExWAS.MLM.bh","Av.ExWAS.MLM.by",
                                                  "Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA"))


ggplot(data = data1.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - Averaged data \nPerformance to identify the true exposure whatever the time point by ICC")


# DLNM

data1.dlnm<-result.data1.100.ICC[result.data1.100.ICC$Method %in% c( "DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                                     "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                                     "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                                     "DLNMselectforward.none","DLNMselectforward.bonf"),]

data1.dlnm$Method<-factor(data1.dlnm$Method, levels=c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                      "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                      "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                      "DLNMselectforward.none","DLNMselectforward.bonf"))


ggplot(data = data1.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - DLM \nPerformance to identify the true exposure whatever the time point by ICC")


# Av + DLNM

data1.dlnm<-result.data1.100.ICC[result.data1.100.ICC$Method %in% c( "Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                                     "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                                     "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by" ),]

data1.dlnm$Method<-factor(data1.dlnm$Method, levels=c("Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                      "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                      "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by"))


ggplot(data = data1.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All time points associated with Y - DLM after Av.ExWAS \nPerformance to identify the true exposure whatever the time point by ICC")




### data2

load(file="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary/summary_data2.100_ICC.2023.Rdata")


result.data2.100.ICC<-result.data2.100.ICC %>% 
  group_by(Method,expo,ICC) %>% 
  summarise(
    false.disc.rate_mean = mean(na.omit(false.disc.rate)),
    sensitivity_mean = mean(na.omit(sensitivity))
  )


result.data2.100.ICC$ICC<-as.character(result.data2.100.ICC$ICC)
result.data2.100.ICC$ICC<-factor(result.data2.100.ICC$ICC, levels=c("1","0.1","0.5","0.9"))

result.data2.100.ICC$color.perf<-as.factor(ifelse(result.data2.100.ICC$sensitivity_mean>80 & result.data2.100.ICC$false.disc.rate_mean<20,"S>80% and FDR<20%",
                                                  ifelse(result.data2.100.ICC$sensitivity_mean>70 & result.data2.100.ICC$false.disc.rate_mean<30,"S>70% and FDR<30%",
                                                         ifelse(result.data2.100.ICC$sensitivity_mean>60 & result.data2.100.ICC$false.disc.rate_mean<40,"S>60% and FDR<40%","S<60% or FDR>40%"))))

# raw data

data2.raw<-result.data2.100.ICC[result.data2.100.ICC$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.MLM.none","Raw.Exwas.bon","Raw.Exwas.MLM.bon","Raw.ExWAS.bh",          
                                                                   "Raw.ExWAS.MLM.bh","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.min","Raw.ENET.opt",          
                                                                   "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"),]

data2.raw$Method<-factor(data2.raw$Method, levels=c("Raw.ExWAS.none","Raw.Exwas.bon","Raw.ExWAS.bh","Raw.ExWAS.by",          
                                                    "Raw.ExWAS.MLM.none","Raw.Exwas.MLM.bon","Raw.ExWAS.MLM.bh","Raw.ExWAS.MLM.by",
                                                    "Raw.ENET.min","Raw.ENET.opt", "Raw.sPLS","Raw.sNPLS","Raw.MMPC","Raw.DSA"))


ggplot(data = data2.raw)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - Raw data \nPerformance to identify the true exposure at the true time point by ICC")



# Av data

data2.av<-result.data2.100.ICC[result.data2.100.ICC$Method %in% c( "Av.ExWAS.none","Av.ExWAS.MLM.none","Av.Exwas.bon",          
                                                                   "Av.Exwas.MLM.bon","Av.ExWAS.bh","Av.ExWAS.MLM.bh","Av.ExWAS.by",           
                                                                   "Av.ExWAS.MLM.by","Av.ENET.min","Av.ENET.opt","Av.sPLS",               
                                                                   "Av.MMPC","Av.DSA"),]

data2.av$Method<-factor(data2.av$Method, levels=c("Av.ExWAS.none","Av.Exwas.bon","Av.ExWAS.bh","Av.ExWAS.by",      
                                                  "Av.ExWAS.MLM.none","Av.Exwas.MLM.bon","Av.ExWAS.MLM.bh","Av.ExWAS.MLM.by",
                                                  "Av.ENET.min","Av.ENET.opt","Av.sPLS","Av.MMPC","Av.DSA"))


ggplot(data = data2.av)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - Averaged data \nPerformance to identify the true exposure at the true time point by ICC")


# DLNM

data2.dlnm<-result.data2.100.ICC[result.data2.100.ICC$Method %in% c( "DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                                     "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                                     "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                                     "DLNMselectforward.none","DLNMselectforward.bonf"),]

data2.dlnm$Method<-factor(data2.dlnm$Method, levels=c("DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                                      "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",
                                                      "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                                      "DLNMselectforward.none","DLNMselectforward.bonf"))


ggplot(data = data2.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - DLM \nPerformance to identify the true exposure at the true time point by ICC")


# Av + DLNM

data2.dlnm<-result.data2.100.ICC[result.data2.100.ICC$Method %in% c( "Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                                     "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                                     "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by" ),]

data2.dlnm$Method<-factor(data2.dlnm$Method, levels=c("Av.DLNMpen.none","Av.DLNMpen.bonf","Av.DLNMpen.bh","Av.DLNMpen.by",
                                                      "Av.DLNMselect.none","Av.DLNMselect.bonf","Av.DLNMselect.bh","Av.DLNMselect.by",
                                                      "Av.DLNMselectback.none","Av.DLNMselectback.bonf","Av.DLNMselectback.bh","Av.DLNMselectback.by"))


ggplot(data = data2.dlnm)+
  geom_point(aes(x=sensitivity_mean,y=false.disc.rate_mean,color=color.perf,shape=ICC),size=2)+
  scale_shape_manual(values=c(16,0,2,5))+
  scale_color_manual(values=c("black","orange","steelblue1","green"))+
  facet_grid(expo~Method)+
  theme_bw()+
  theme(strip.text = element_text(size = 8, margin = margin()),legend.position = "bottom",legend.direction = "horizontal",
        axis.text = element_text(size=6))+
  xlim(0,100)+
  ylim(0,100)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y - DLM after Av.ExWAS \nPerformance to identify the true exposure at the true time point by ICC")

