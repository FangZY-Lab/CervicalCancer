####### GSE7803
############## 1 limma #########################################################
############## Loading Library and data ########################################

rm(list = ls())
gc()
library(limma)
setwd("/Users/wuhongbo/Desktop/Cervical Cancer/data")
g=read.csv("GSE7803_group.csv",header=T,row.names=1)
Data=read.csv("GSE7803_Data.csv",header=T,row.names=1)
dim(Data)
############ linear model using limma ##########################################
groupid <- ifelse(g$Group%in%c('control'), 1,2)
design <- model.matrix(~ -1+factor(groupid))
colnames(design) <- c("Control", "Case")  
limma.fit <- lmFit(log(Data,2), design)
contrast.matrix <- makeContrasts(Control-Case, levels=design)
limma.fit2 <- contrasts.fit(limma.fit, contrast.matrix)
limma.fit2 <- eBayes(limma.fit2)
Pvalue_limma <- p.adjust(limma.fit2$p.value, method="BH")
######### Differential expressed genes using adj.p value  ######################
data.adj.p=data.frame(limma.fit2$p.value,Pvalue_limma)
colnames(data.adj.p) <- c("p_val", "adj_pval")
sig.genes=data.adj.p[(which(data.adj.p$adj_pval< 0.05)),]
sig.genesordered=sig.genes[order(sig.genes$adj_pval),]
dim(sig.genes)
write.csv(as.data.frame(data.adj.p), file="limma_GSE7803.csv")
write.csv(as.data.frame(sig.genesordered), file="limma_GSE7803_DEG.csv")

################################################################################
################## 2  SAM Analysis  ############################################

rm(list = ls())
gc()
library(samr)
g=read.csv("GSE7803_group.csv",header=T,row.names=1)
Data0=as.matrix(read.csv("GSE7803_Data.csv",header=T,row.names=1))
head(Data0)
###  process the data for samr #################################################
y1 <- ifelse(g$Group%in%c('control'), 1,2)
data=list(x=Data0,y=y1,logged2=TRUE)
samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)
pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
P_adj <- p.adjust(pv, method="BH")
#### Differentially Expressed Genes ####################################
data.adj.p=data.frame(samr.obj$tt,samr.obj$ttstar,P_adj)
sig.genes=data.adj.p[(which(data.adj.p$P_adj< 0.05)),]
sig.genesordered=sig.genes[order(sig.genes$P_adj),102,drop=F]
write.csv(as.data.frame(data.adj.p), file="SAM_GSE7803.csv")
write.csv(as.data.frame(sig.genesordered), file="SAM_GSE7803_DEG.csv")

################################################################################
#######  3 t/FCROS Analysis   ########################################

rm(list = ls())
gc()
library(fcros)
Data=read.csv("GSE7803_Data.csv",header=T,row.names=1)
dim(Data)
cont=c("GSM189384","GSM189385","GSM189386","GSM189387","GSM189388","GSM189389","GSM189390","GSM189391","GSM189392","GSM189393")
test=c("GSM189381","GSM189382","GSM189383","GSM189394","GSM189395","GSM189396","GSM189397","GSM189398","GSM189399","GSM189400","GSM189401","GSM189402","GSM189403","GSM189404","GSM189405","GSM189406","GSM189407","GSM189408","GSM189409","GSM189410","GSM189411","GSM189412","GSM189413","GSM189414","GSM189415","GSM189416","GSM189417","GSM189418","GSM189419","GSM189420","GSM189421")
log2.opt <- 0
at <- fcrosTtest(Data, cont, test, log2.opt)
#####  Differentially expressed genes by using p adjusted value ##################
adj_pval <- p.adjust(at$p.value, method="BH")
dt.f=data.frame(at$idnames,adj_pval)
sig.genes=dt.f[(which(dt.f$adj_pval< 0.05)),]
dim(sig.genes)
sig.genesordered=sig.genes[order(sig.genes$adj_pval),]
rownames(sig.genes)=sig.genes$at.idnames
rownames(sig.genesordered)=sig.genesordered$at.idnames
sig.genes=sig.genes[,-1,drop=F]
sig.genesordered=sig.genesordered[,-1,drop=F]
write.csv(as.data.frame(dt.f), file="t_GSE7803.csv")
write.csv(as.data.frame(sig.genesordered), file="t_GSE7803_DEG.csv")



######## GSE9750
############## 1 limma #########################################################
############## Loading Library and data ########################################

rm(list = ls())
gc()
library(limma)
g=read.csv("GSE9750_group.csv",header=T,row.names=1)
Data=read.csv("GSE9750_Data.csv",header=T,row.names=1)
dim(Data)
############ linear model using limma ##########################################
groupid <- ifelse(g$Group%in%c('control'), 1,2)
design <- model.matrix(~ -1+factor(groupid))
colnames(design) <- c("Control", "Case")  
limma.fit <- lmFit(log(Data,2), design)
contrast.matrix <- makeContrasts(Control-Case, levels=design)
limma.fit2 <- contrasts.fit(limma.fit, contrast.matrix)
limma.fit2 <- eBayes(limma.fit2)
Pvalue_limma <- p.adjust(limma.fit2$p.value, method="BH")
######### Differential expressed genes using adj.p value  ######################
data.adj.p=data.frame(limma.fit2$p.value,Pvalue_limma)
colnames(data.adj.p) <- c("p_val", "adj_pval")
sig.genes=data.adj.p[(which(data.adj.p$adj_pval< 0.05)),]
sig.genesordered=sig.genes[order(sig.genes$adj_pval),]
dim(sig.genes)
write.csv(as.data.frame(data.adj.p), file="limma_GSE9750.csv")
write.csv(as.data.frame(sig.genesordered), file="limma_GSE9750_DEG.csv")

################################################################################
################## 2  SAM Analysis  ############################################

rm(list = ls())
gc()
library(samr)
g=read.csv("GSE9750_group.csv",header=T,row.names=1)
Data0=as.matrix(read.csv("GSE9750_Data.csv",header=T,row.names=1))
head(Data0)
###  process the data for samr #################################################
y1 <- ifelse(g$Group%in%c('control'), 1,2)
data=list(x=Data0,y=y1,logged2=TRUE)
samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)
pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
P_adj <- p.adjust(pv, method="BH")
#### Differentially Expressed Genes ####################################
data.adj.p=data.frame(samr.obj$tt,samr.obj$ttstar,P_adj)
#colnames(data.adj.p) <- c("p_val", "adj_pval")
sig.genes=data.adj.p[(which(data.adj.p$P_adj< 0.05)),]
sig.genesordered=sig.genes[order(sig.genes$P_adj),102,drop=F]
write.csv(as.data.frame(data.adj.p), file="SAM_GSE9750.csv")
write.csv(as.data.frame(sig.genesordered), file="SAM_GSE9750_DEG.csv")

################################################################################
#######  3 t/FCROS Analysis   ########################################

rm(list = ls())
gc()
library(fcros)
g=read.csv("GSE9750_group.csv",header=T,row.names=1)
Data=read.csv("GSE9750_Data.csv",header=T,row.names=1)
dim(Data)
cont=rownames(g[g$Group == "control", ])
test=rownames(g[g$Group == "case", ])
log2.opt <- 0
at <- fcrosTtest(Data, cont, test, log2.opt)
#####  Differentially expressed genes by using p adjusted value ##################
adj_pval <- p.adjust(at$p.value, method="BH")
dt.f=data.frame(at$idnames,adj_pval)
sig.genes=dt.f[(which(dt.f$adj_pval< 0.05)),]
dim(sig.genes)
sig.genesordered=sig.genes[order(sig.genes$adj_pval),]
rownames(sig.genes)=sig.genes$at.idnames
rownames(sig.genesordered)=sig.genesordered$at.idnames
sig.genes=sig.genes[,-1,drop=F]
sig.genesordered=sig.genesordered[,-1,drop=F]
write.csv(as.data.frame(dt.f), file="t_GSE9750.csv")
write.csv(as.data.frame(sig.genesordered), file="t_GSE9750_DEG.csv")

