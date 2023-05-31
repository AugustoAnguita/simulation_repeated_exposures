library(ggplot2)

setwd(dir="D:/Home/cwarembourg/Documents/Lifecycle/simulation mars 2022/Summary")


load("summary.data1.2023.Rdata")
load("summary.data2.2023.Rdata")
load("summary.data1.100.2023.Rdata")
load("summary.data2.100.2023.Rdata")


unique(summary_data1$Method)

### data1 at true time point

data1<-summary_data1[summary_data1$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.DSA",
                                                 "Av.ExWAS.none","Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.opt","Av.sPLS","Av.DSA",                  
                                                 "DLNMselect.by","DLNMselect.AVG.by"),]

data1$type_data<-ifelse(data1$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.DSA","DLNMselect.by"),"Raw","Averaged")

data1$false.disc.rate_mean_3<-as.numeric(data1$false.disc.rate_mean_3)
data1$false.disc.rate_mean_5<-as.numeric(data1$false.disc.rate_mean_5)
data1$false.disc.rate_mean_10<-as.numeric(data1$false.disc.rate_mean_10)
data1$sensitivity_mean_3<-as.numeric(data1$sensitivity_mean_3)
data1$sensitivity_mean_5<-as.numeric(data1$sensitivity_mean_5)
data1$sensitivity_mean_10<-as.numeric(data1$sensitivity_mean_10)


levels(data1$Method) <- list(ExWAS.none="Raw.ExWAS.none", ExWAS.by="Raw.ExWAS.by", ExWAS.MLM.by="Raw.ExWAS.MLM.by", ENET="Raw.ENET.opt", 
                             sPLS="Raw.sPLS", sNPLS="Raw.sNPLS", DSA="Raw.DSA", DLNMselect.by="DLNMselect.by", 
                             ExWAS.none="Av.ExWAS.none", ExWAS.by="Av.ExWAS.by", ExWAS.MLM.by="Av.ExWAS.MLM.by", ENET="Av.ENET.opt", 
                             sPLS="Av.sPLS", DSA="Av.DSA", DLNMselect.by="DLNMselect.AVG.by")

data1$type_data <- factor(data1$type_data, levels = c("Raw","Averaged"))


## n=3
ggplot(data=data1,aes(x=sensitivity_mean_3,y=false.disc.rate_mean_3,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+

  xlim(0,100)+
  ylim(0,100)+
  # scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100))+
  # scale_y_continuous(breaks=c(0, 20, 40, 60, 80, 100))+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All 5 time points are associated with Y \n n=3 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=5
ggplot(data=data1,aes(x=sensitivity_mean_5,y=false.disc.rate_mean_5,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All 5 time points are associated with Y \n n=5 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=10
ggplot(data=data1,aes(x=sensitivity_mean_10,y=false.disc.rate_mean_10,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All 5 time points are associated with Y \n n=10 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")






### data2 at true time point

data2<-summary_data2[summary_data2$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.opt","Raw.sPLS","Raw.sNPLS","Raw.DSA",
                                                 "Av.ExWAS.none","Av.ExWAS.by","Av.ExWAS.MLM.by","Av.ENET.opt","Av.sPLS","Av.DSA",                  
                                                 "DLNMselect.by","DLNMselect.AVG.by"),]

data2$type_data<-ifelse(data2$Method %in% c("Raw.ExWAS.none","Raw.ExWAS.by","Raw.ExWAS.MLM.by","Raw.ENET.opt","Raw.sPLS","Raw.nSPLS","Raw.DSA","DLNMselect.by"),"Raw","Averaged")

data2$false.disc.rate_mean_3<-as.numeric(data2$false.disc.rate_mean_3)
data2$false.disc.rate_mean_5<-as.numeric(data2$false.disc.rate_mean_5)
data2$false.disc.rate_mean_10<-as.numeric(data2$false.disc.rate_mean_10)
data2$sensitivity_mean_3<-as.numeric(data2$sensitivity_mean_3)
data2$sensitivity_mean_5<-as.numeric(data2$sensitivity_mean_5)
data2$sensitivity_mean_10<-as.numeric(data2$sensitivity_mean_10)


levels(data2$Method) <- list(ExWAS.none="Raw.ExWAS.none", ExWAS.by="Raw.ExWAS.by", ExWAS.MLM.by="Raw.ExWAS.MLM.by", ENET="Raw.ENET.opt", 
                             sPLS="Raw.sPLS",sNPLS="Raw.sNPLS", DSA="Raw.DSA", DLNMselect.by="DLNMselect.by", 
                             ExWAS.none="Av.ExWAS.none", ExWAS.by="Av.ExWAS.by", ExWAS.MLM.by="Av.ExWAS.MLM.by", ENET="Av.ENET.opt", 
                             sPLS="Av.sPLS", DSA="Av.DSA", DLNMselect.by="DLNMselect.AVG.by")

data2$type_data <- factor(data2$type_data, levels = c("Raw","Averaged"))


## n=3
ggplot(data=data2,aes(x=sensitivity_mean_3,y=false.disc.rate_mean_3,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y \n n=3 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=5
ggplot(data=data2,aes(x=sensitivity_mean_5,y=false.disc.rate_mean_5,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y \n n=5 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=10
ggplot(data=data2,aes(x=sensitivity_mean_10,y=false.disc.rate_mean_10,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y \n n=10 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")






### data1 whatever true time point

data1.100<-summary_data1.100[summary_data1.100$Method %in% c("Raw.EWAS.none","Raw.EWAS.by","Raw.EWAS_LM.by","Raw.ENET_OPT","Raw.SPLS","Raw.sNPLS","Raw.DSA",
                                                 "Av.2step.EWAS.none","Av.2step.EWAS.by","Av.2step.EWAS_LM.by","Av.2step.ENET_OPT","Av.2step.SPLS","Av.2step.DSA",                  
                                                 "DLNMselect.2step.by","DLNMselect.AVG.2step.by"),]
 
data1.100$type_data<-ifelse(data1.100$Method %in% c("Raw.EWAS.none","Raw.EWAS.by","Raw.EWAS.LM.by","Raw.ENET_OPT","Raw.SPLS","Raw.sNPLS","Raw.DSA","DLNMselect.2step.by"),"Raw","Averaged")

data1.100$false.disc.rate_mean_3<-as.numeric(data1.100$false.disc.rate_mean_3)
data1.100$false.disc.rate_mean_5<-as.numeric(data1.100$false.disc.rate_mean_5)
data1.100$false.disc.rate_mean_10<-as.numeric(data1.100$false.disc.rate_mean_10)
data1.100$sensitivity_mean_3<-as.numeric(data1.100$sensitivity_mean_3)
data1.100$sensitivity_mean_5<-as.numeric(data1.100$sensitivity_mean_5)
data1.100$sensitivity_mean_10<-as.numeric(data1.100$sensitivity_mean_10)


levels(data1.100$Method) <- list(ExWAS.none="Raw.EWAS.none", ExWAS.by="Raw.EWAS.by", ExWAS.MLM.by="Raw.EWAS_LM.by", ENET="Raw.ENET_OPT", 
                             sPLS="Raw.SPLS",sNPLS="Raw.sNPLS", DSA="Raw.DSA", DLNMselect.by="DLNMselect.2step.by", 
                             ExWAS.none="Av.2step.EWAS.none", ExWAS.by="Av.2step.EWAS.by", ExWAS.MLM.by="Av.2step.EWAS_LM.by", ENET="Av.2step.ENET_OPT", 
                             sPLS="Av.2step.SPLS", DSA="Av.2step.DSA", DLNMselect.by="DLNMselect.AVG.2step.by")

data1.100$type_data <- factor(data1.100$type_data, levels = c("Raw","Averaged"))


## n=3
ggplot(data=data1.100,aes(x=sensitivity_mean_3,y=false.disc.rate_mean_3,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All 5 time points are associated with Y \n n=3 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=5
ggplot(data=data1.100,aes(x=sensitivity_mean_5,y=false.disc.rate_mean_5,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All 5 time points are associated with Y \n n=5 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=10
ggplot(data=data1.100,aes(x=sensitivity_mean_10,y=false.disc.rate_mean_10,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 1: All 5 time points are associated with Y \n n=10 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


### data2 whatever true time point

data2.100<-summary_data2.100[summary_data2.100$Method %in% c("Raw.EWAS.none","Raw.EWAS.by","Raw.EWAS_LM.by","Raw.ENET_OPT","Raw.SPLS","Raw.sNPLS","Raw.DSA",
                                                             "Av.2step.EWAS.none","Av.2step.EWAS.by","Av.2step.EWAS_LM.by","Av.2step.ENET_OPT","Av.2step.SPLS","Av.2step.DSA",                  
                                                             "DLNMselect.2step.by","DLNMselect.AVG.2step.by"),]

data2.100$type_data<-ifelse(data2.100$Method %in% c("Raw.EWAS.none","Raw.EWAS.by","Raw.EWAS.LM.by","Raw.ENET_OPT","Raw.SPLS","Raw.sNPLS","Raw.DSA","DLNMselect.2step.by"),"Raw","Averaged")

data2.100$false.disc.rate_mean_3<-as.numeric(data2.100$false.disc.rate_mean_3)
data2.100$false.disc.rate_mean_5<-as.numeric(data2.100$false.disc.rate_mean_5)
data2.100$false.disc.rate_mean_10<-as.numeric(data2.100$false.disc.rate_mean_10)
data2.100$sensitivity_mean_3<-as.numeric(data2.100$sensitivity_mean_3)
data2.100$sensitivity_mean_5<-as.numeric(data2.100$sensitivity_mean_5)
data2.100$sensitivity_mean_10<-as.numeric(data2.100$sensitivity_mean_10)


levels(data2.100$Method) <- list(ExWAS.none="Raw.EWAS.none", ExWAS.by="Raw.EWAS.by", ExWAS.MLM.by="Raw.EWAS_LM.by", ENET="Raw.ENET_OPT", 
                                 sPLS="Raw.SPLS",sNPLS="Raw.sNPLS", DSA="Raw.DSA", DLNMselect.by="DLNMselect.2step.by", 
                                 ExWAS.none="Av.2step.EWAS.none", ExWAS.by="Av.2step.EWAS.by", ExWAS.MLM.by="Av.2step.EWAS_LM.by", ENET="Av.2step.ENET_OPT", 
                                 sPLS="Av.2step.SPLS", DSA="Av.2step.DSA", DLNMselect.by="DLNMselect.AVG.2step.by")

data2.100$type_data <- factor(data2.100$type_data, levels = c("Raw","Averaged"))


## n=3
ggplot(data=data2.100,aes(x=sensitivity_mean_3,y=false.disc.rate_mean_3,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y \n n=3 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=5
ggplot(data=data2.100,aes(x=sensitivity_mean_5,y=false.disc.rate_mean_5,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y \n n=5 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")


## n=10
ggplot(data=data2.100,aes(x=sensitivity_mean_10,y=false.disc.rate_mean_10,col=type_data,group=Method,shape=Method))+
  geom_point(size=4)+
  scale_shape_manual(values=c(15,17,18,19,3,4,8,9))+
  scale_color_manual(values=c("black","red"))+
  
  xlim(0,100)+
  ylim(0,100)+
  geom_text(aes(label=Method),nudge_y=3,check_overlap = T)+
  xlab("Sensitivity")+
  ylab("False Discovery Rate")+
  ggtitle("Scenario 2: A single time point is associated with Y \n n=10 true exposures")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size=20), 
        legend.position = "bottom",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=16))+
  guides(shape = "none")



