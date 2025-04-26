####### GSE168244
############## 1 DEseq2 ########################################################
############## Loading Library and data ########################################

rm(list = ls())
gc()
library(apeglm)
library(DESeq2)
library(readxl)
dat=read.csv("~/Desktop/Cervical Cancer/data/GSE168244_data.csv",header=TRUE,row.names=1)
coldata=read.csv("~/Desktop/Cervical Cancer/coldata/GSE168244_coldata.csv",header=TRUE,sep = ',')
dim(dat)
cts=as.matrix(dat)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
######################## DEseq2 Analysis #######################################
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >0
dds <- dds[keep,]
ddsDE <- DESeq(dds)
resultsNames(ddsDE)
res <- results(ddsDE, name="condition_Control_vs_Case")
summary(res)
res05 <- res[which(res$padj<0.05), ]
dim(res05)
res05Ordered <- res05[order(res05$padj),]
write.csv(as.data.frame(res), file="DESeq2_GSE168244.csv")
write.csv(as.data.frame(res05Ordered), file="DESeq2_GSE168244_DEG.csv")

############## 2 edgeR #########################################################
############## Loading Library and data ########################################
rm(list = ls())
gc()
library(edgeR)
Data=read.csv("~/Desktop/Cervical Cancer/data/GSE168244_data.csv",header=TRUE,row.names=1)
coldata=read.csv("~/Desktop/Cervical Cancer/coldata/GSE168244_coldata.csv",header=TRUE,sep = ',')
dim(Data)
######################## edgeR Analysis #######################################
n1=18
n2=17
groupid<-c(rep(1,n1),rep(2,n2))
y=DGEList(counts=Data,group=groupid)
keep <- rowSums(cpm(y) > 0) >= 2
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ 0+factor(groupid))
colnames(design) <- c("Control", "Case")  
library(statmod)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
summary(fit$df.prior)
B.LvsP=makeContrasts(Control-Case, levels=design)
res <- glmQLFTest(fit, contrast=B.LvsP)
adj_pval <- p.adjust(res$table$PValue, method="BH")
data.adj.p=data.frame(res$table,adj_pval)
sig.genes=data.adj.p[(which(data.adj.p$adj_pval< 0.05)),]
sig.genesordered=sig.genes[order(sig.genes$adj_pval),]
dim(sig.genes)
write.csv(as.data.frame(data.adj.p), file="edgeR_GSE168244.csv")
write.csv(as.data.frame(sig.genesordered), file="edgeR_GSE168244_DEG.csv")

############## 3 NOISeq #########################################################
############## Loading Library and data ########################################
rm(list = ls())
gc()
library(NOISeq)
Data=read.csv("~/Desktop/Cervical Cancer/data/GSE168244_data.csv",header=TRUE,row.names=1)
coldata=read.csv("~/Desktop/Cervical Cancer/coldata/GSE168244_coldata.csv",header=TRUE,sep = ',')
dat=Data[rowSums(Data[, -1])>2,]
dim(dat)
dim(coldata)
############### NOISeq Analysis #################################################
mydata=readData(data = dat,factors = coldata)
noi=noiseqbio(mydata, k = 0.5, norm = "tmm", factor="condition",  adj = 0.05,lc = 0,
              a0per = 0.9,conditions=c("Control","Case"))
mynoiseq.deg = degenes(noi, q = 0.8, M = NULL)
write.csv(as.data.frame(mynoiseq.deg), file="NOISeq_GSE168244_DEG.csv")


####### GSE149450
############## 1 DEseq2 ########################################################
############## Loading Library and data ########################################

rm(list = ls())
gc()
library(apeglm)
library(DESeq2)
library(readxl)
dat=read.csv("~/Desktop/Cervical Cancer/data/GSE149450_data.csv",header=TRUE,row.names=1)
coldata=read.csv("~/Desktop/Cervical Cancer/coldata/GSE149450_coldata.csv",header=TRUE,sep = ',')
dim(dat)
cts=as.matrix(dat)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
######################## DEseq2 Analysis #######################################
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >0
dds <- dds[keep,]
ddsDE <- DESeq(dds)
resultsNames(ddsDE)
res <- results(ddsDE, name="condition_Control_vs_Case")
summary(res)
res05 <- res[which(res$padj<0.05), ]
dim(res05)
res05Ordered <- res05[order(res05$padj),]
write.csv(as.data.frame(res), file="DESeq2_GSE149450.csv")
write.csv(as.data.frame(res05Ordered), file="DESeq2_GSE149450_DEG.csv")

############## 2 edgeR #########################################################
############## Loading Library and data ########################################
rm(list = ls())
gc()
library(edgeR)
Data=read.csv("~/Desktop/Cervical Cancer/data/GSE149450_data.csv",header=TRUE,row.names=1)
dim(Data)
######################## edgeR Analysis #######################################
n1=2
n2=2
groupid<-c(rep(1,n1),rep(2,n2))
y=DGEList(counts=Data,group=groupid)
keep <- rowSums(cpm(y) > 0) >= 2
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ 0+factor(groupid))
colnames(design) <- c("Control", "Case")  
library(statmod)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
summary(fit$df.prior)
B.LvsP=makeContrasts(Control-Case, levels=design)
res <- glmQLFTest(fit, contrast=B.LvsP)
adj_pval <- p.adjust(res$table$PValue, method="BH")
data.adj.p=data.frame(res$table,adj_pval)
sig.genes=data.adj.p[(which(data.adj.p$adj_pval< 0.05)),]
sig.genesordered=sig.genes[order(sig.genes$adj_pval),]
dim(sig.genes)
write.csv(as.data.frame(data.adj.p), file="edgeR_GSE149450.csv")
write.csv(as.data.frame(sig.genesordered), file="edgeR_GSE149450_DEG.csv")

############## 3 NOISeq #########################################################
############## Loading Library and data ########################################
rm(list = ls())
gc()
library(lfc)
library(NOISeq)
Data=read.csv("~/Desktop/Cervical Cancer/data/GSE149450_data.csv",header=TRUE,row.names=1)
coldata=read.csv("~/Desktop/Cervical Cancer/coldata/GSE149450_coldata.csv",header=TRUE,sep = ',')
dat=Data[rowSums(Data[, -1])>2,]
dim(dat)
dim(coldata)
############### NOISeq Analysis #################################################
mydata=readData(data = dat,factors = coldata)
noi=noiseqbio(mydata, k = 0.5, norm = "tmm", factor="condition",lc = 0, r = 50, plot = FALSE,
              a0per = 0.9,filter=1,conditions=c("Control","Case"),cv.cutoff = 100,adj = 0.05,cpm = 1)
mynoiseq.deg = degenes(noi, q = 0.965, M = NULL)
write.csv(as.data.frame(mynoiseq.deg), file="NOISeq_GSE149450_DEG.csv")

