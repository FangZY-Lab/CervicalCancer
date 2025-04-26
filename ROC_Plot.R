setwd("G:\\")

library(pROC)
#=====================Data_GSE7803=======================
pdf("ROC_GSE7803.pdf")
Data<-read.csv("ROCdata_GSE7803.csv", header=T)
n1=10; n2=31
rownames(Data)<-Data[,1];Data<-Data[,-1];Data; dim(Data)
Data<-log2(Data)
df <- data.frame(t(Data)); labels<-c(rep(0,n1), rep(1,n2))
df$labels<-labels
#======================================================
#========PLOT=========================================
#par(mfrow = c(2, 2)) 
par(pty="s")
roc1 <- roc(labels ~ CDK1, df, smooth=TRUE, print.auc = TRUE)
plot(roc1, show.thres=TRUE, legacy.axes=TRUE,col="#1f77b4",lwd=4)
auc_roc1<- as.vector(auc(roc1))            

roc2 <- roc(labels ~ TOP2A, df, smooth=TRUE,print.auc = TRUE)
plot(roc2, add=TRUE, percent=roc1$percent, col="#2ca02c",lwd=4)
auc_roc2<- as.vector(auc(roc2)) 

roc3<- roc(labels ~ BRCA1, df, smooth=TRUE, print.auc = TRUE)
plot(roc3, add=TRUE, percent=roc1$percent, col="#008080",lwd=4)
auc_roc3<- as.vector(auc(roc3)) 

roc4<- roc(labels ~ CDC20, df, smooth=TRUE,print.auc = TRUE)
plot(roc4, add=TRUE, percent=roc1$percent, col="#ff7f0e",lwd=4)
auc_roc4<- as.vector(auc(roc4)) 


roc5 <- roc(labels ~ FEN1, df, smooth=TRUE,print.auc = TRUE)
plot(roc5, add=TRUE, percent=roc1$percent, col="#9467bd",lwd=4)
auc_roc5<- as.vector(auc(roc5)) 


roc6<- roc(labels ~ BUB1B, df, smooth=TRUE, print.auc = TRUE)
plot(roc6, add=TRUE, percent=roc1$percent, col="#bcbd22",lwd=4)
auc_roc6<- as.vector(auc(roc6)) 

roc7<- roc(labels ~ CDC45, df, smooth=TRUE,print.auc = TRUE)
plot(roc7, add=TRUE, percent=roc1$percent, col="#e377c2",lwd=4)
auc_roc7<- as.vector(auc(roc7)) 

roc8<- roc(labels ~ KIF23, df, smooth=TRUE, print.auc = TRUE)
plot(roc8, add=TRUE, percent=roc1$percent, col="#17becf",lwd=4)
auc_roc8<- as.vector(auc(roc8)) 


roc9 <- roc(labels ~ KIF11, df, smooth=TRUE,print.auc = TRUE)
plot(roc9,add=TRUE, percent=roc1$percent, col="#8c564b",lwd=4)
auc_roc9<- as.vector(auc(roc9)) 

roc10<- roc(labels ~ MCM7, df, smooth=TRUE, print.auc = TRUE)
plot(roc10, add=TRUE, percent=roc1$percent, col="#ffbb78",lwd=4)
auc_roc10<- as.vector(auc(roc10)) 

#=====================AUC
AUC<-round(mean(c(auc_roc1,auc_roc2,auc_roc3,auc_roc4,
                  auc_roc5,auc_roc6,auc_roc7,auc_roc8,auc_roc9,auc_roc10)),digit=4)
#=========legend
#legend("bottomright", legend=c("S100B", "NDKA"), col=c("#1c61b6", "#008600"), lwd=2)
text(.28, .15, labels=paste("Avg. AUC =", AUC), adj=c(0,.5), col = "1",font = 2)
text(.28, .10, labels=paste("GEO ID: GSE7803"), adj=c(0,.5),col = "1",font = 2)

dev.off()


#********************************************************************
#********************************************************************
#********************************************************************
#=======================================================
library(pROC)
#=====================Data_GSE9750=======================
pdf("ROC_GSE9750.pdf")
Data<-read.csv("ROCdata_GSE9750.csv", header=T)
n1=24; n2=33
rownames(Data)<-Data[,1];Data<-Data[,-c(1:10)];Data; dim(Data)
Data<-log2(Data)
df <- data.frame(t(Data)); labels<-c(rep(0,n1), rep(1,n2))
df$labels<-labels
#======================================================
#========PLOT=========================================
#par(mfrow = c(2, 2)) 
par(pty="s")
roc1 <- roc(labels ~ CDK1, df, smooth=TRUE, print.auc = TRUE)
plot(roc1, show.thres=TRUE, legacy.axes=TRUE,col="#1f77b4",lwd=4)
auc_roc1<- as.vector(auc(roc1))            

roc2 <- roc(labels ~ TOP2A, df, smooth=TRUE,print.auc = TRUE)
plot(roc2, add=TRUE, percent=roc1$percent, col="#2ca02c",lwd=4)
auc_roc2<- as.vector(auc(roc2)) 

roc3<- roc(labels ~ BRCA1, df, smooth=TRUE, print.auc = TRUE)
plot(roc3, add=TRUE, percent=roc1$percent, col="#008080",lwd=4)
auc_roc3<- as.vector(auc(roc3)) 

roc4<- roc(labels ~ CDC20, df, smooth=TRUE,print.auc = TRUE)
plot(roc4, add=TRUE, percent=roc1$percent, col="#ff7f0e",lwd=4)
auc_roc4<- as.vector(auc(roc4)) 


roc5 <- roc(labels ~ FEN1, df, smooth=TRUE,print.auc = TRUE)
plot(roc5, add=TRUE, percent=roc1$percent, col="#9467bd",lwd=4)
auc_roc5<- as.vector(auc(roc5)) 


roc6<- roc(labels ~ BUB1B, df, smooth=TRUE, print.auc = TRUE)
plot(roc6, add=TRUE, percent=roc1$percent, col="#bcbd22",lwd=4)
auc_roc6<- as.vector(auc(roc6)) 

roc7<- roc(labels ~ CDC45, df, smooth=TRUE,print.auc = TRUE)
plot(roc7, add=TRUE, percent=roc1$percent, col="#e377c2",lwd=4)
auc_roc7<- as.vector(auc(roc7)) 

roc8<- roc(labels ~ KIF23, df, smooth=TRUE, print.auc = TRUE)
plot(roc8, add=TRUE, percent=roc1$percent, col="#17becf",lwd=4)
auc_roc8<- as.vector(auc(roc8)) 


roc9 <- roc(labels ~ KIF11, df, smooth=TRUE,print.auc = TRUE)
plot(roc9,add=TRUE, percent=roc1$percent, col="#8c564b",lwd=4)
auc_roc9<- as.vector(auc(roc9)) 

roc10<- roc(labels ~ MCM7, df, smooth=TRUE, print.auc = TRUE)
plot(roc10, add=TRUE, percent=roc1$percent, col="#ffbb78",lwd=4)
auc_roc10<- as.vector(auc(roc10)) 

#=====================AUC
AUC<-round(mean(c(auc_roc1,auc_roc2,auc_roc3,auc_roc4,
                  auc_roc5,auc_roc6,auc_roc7,auc_roc8,auc_roc9,auc_roc10)),digit=4)
#=========legend
#legend("bottomright", legend=c("S100B", "NDKA"), col=c("#1c61b6", "#008600"), lwd=2)
text(.28, .15, labels=paste("Avg. AUC =", AUC), adj=c(0,.5), col = "1",font = 2)
text(.28, .10, labels=paste("GEO ID: GSE9750"), adj=c(0,.5),col = "1",font = 2)

dev.off()


#********************************************************************
#********************************************************************
#********************************************************************
#=======================================================
library(pROC)
#=====================Data_GSE168244=======================
pdf("ROC_GSE168244.pdf")
Data<-read.csv("ROCdata_GSE168244.csv", header=T)
n1=18; n2=17
rownames(Data)<-Data[,1];Data<-Data[,-1];Data; dim(Data)
#Data<-log2(Data)
df <- data.frame(t(Data)); labels<-c(rep(0,n1), rep(1,n2))
df$labels<-labels
#======================================================
#========PLOT=========================================
#par(mfrow = c(2, 2)) 
par(pty="s")
roc1 <- roc(labels ~ CDK1, df, smooth=TRUE, print.auc = TRUE)
plot(roc1, show.thres=TRUE, legacy.axes=TRUE,col="#1f77b4",lwd=4)
auc_roc1<- as.vector(auc(roc1))            

roc2 <- roc(labels ~ TOP2A, df, smooth=TRUE,print.auc = TRUE)
plot(roc2, add=TRUE, percent=roc1$percent, col="#2ca02c",lwd=4)
auc_roc2<- as.vector(auc(roc2)) 

roc3<- roc(labels ~ BRCA1, df, smooth=TRUE, print.auc = TRUE)
plot(roc3, add=TRUE, percent=roc1$percent, col="#008080",lwd=4)
auc_roc3<- as.vector(auc(roc3)) 

roc4<- roc(labels ~ CDC20, df, smooth=TRUE,print.auc = TRUE)
plot(roc4, add=TRUE, percent=roc1$percent, col="#ff7f0e",lwd=4)
auc_roc4<- as.vector(auc(roc4)) 


roc5 <- roc(labels ~ FEN1, df, smooth=TRUE,print.auc = TRUE)
plot(roc5, add=TRUE, percent=roc1$percent, col="#9467bd",lwd=4)
auc_roc5<- as.vector(auc(roc5)) 


roc6<- roc(labels ~ BUB1B, df, smooth=TRUE, print.auc = TRUE)
plot(roc6, add=TRUE, percent=roc1$percent, col="#bcbd22",lwd=4)
auc_roc6<- as.vector(auc(roc6)) 

roc7<- roc(labels ~ CDC45, df, smooth=TRUE,print.auc = TRUE)
plot(roc7, add=TRUE, percent=roc1$percent, col="#e377c2",lwd=4)
auc_roc7<- as.vector(auc(roc7)) 

roc8<- roc(labels ~ KIF23, df, smooth=TRUE, print.auc = TRUE)
plot(roc8, add=TRUE, percent=roc1$percent, col="#17becf",lwd=4)
auc_roc8<- as.vector(auc(roc8)) 


roc9 <- roc(labels ~ KIF11, df, smooth=TRUE,print.auc = TRUE)
plot(roc9,add=TRUE, percent=roc1$percent, col="#8c564b",lwd=4)
auc_roc9<- as.vector(auc(roc9)) 

roc10<- roc(labels ~ MCM7, df, smooth=TRUE, print.auc = TRUE)
plot(roc10, add=TRUE, percent=roc1$percent, col="#ffbb78",lwd=4)
auc_roc10<- as.vector(auc(roc10)) 

#=====================AUC
AUC<-round(mean(c(auc_roc1,auc_roc2,auc_roc3,auc_roc4,
                  auc_roc5,auc_roc6,auc_roc7,auc_roc8,auc_roc9,auc_roc10)),digit=4)
#=========legend
#legend("bottomright", legend=c("S100B", "NDKA"), col=c("#1c61b6", "#008600"), lwd=2)
text(.28, .15, labels=paste("Avg. AUC =", AUC), adj=c(0,.5), col = "1",font = 2)
text(.28, .10, labels=paste("GEO ID: GSE168244"), adj=c(0,.5),col = "1",font = 2)

dev.off()


