######GSE7803
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*600)
gse1<-getGEO("GSE7803",destdir = ".",getGPL = T,AnnotGPL = T)
cli<-pData(gse1[[1]])
Data=read.csv("GSE7803_Data.csv",header=T,row.names=1)
g=cli[colnames(Data),c(2,10)]
colnames(g)=c('Group','Characteristics')
g$Group=ifelse(g$Characteristics%in%c('normal'),'control','case')
write.csv(g, file="GSE7803_group.csv")

######GSE9750
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*600)
gse1<-getGEO("GSE9750",destdir = ".",getGPL = T,AnnotGPL = T)
cli<-pData(gse1[[1]])
Data=read.csv("GSE9750_Data.csv",header=T,row.names=1)
g=cli[colnames(Data),c(1,8)]
colnames(g)=c('Group','Characteristics')
g$Characteristics=sub("^(.*?)\\s.*$", "\\1", g$Group)
g$Group=ifelse(g$Characteristics%in%c('Normal'),'control','case')
write.csv(g, file="GSE9750_group.csv")

