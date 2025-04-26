############################################################ GSE7803 fgsea ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"
#gmtdir='/Users/wuhongbo/Desktop/数据分析/IPF/fgsea/output'

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                             category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
limma_GSE7803_DEG = read.table('limma_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE7803_DEG = read.table('SAM_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE7803_DEG = read.table('t_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE7803_DEG <- limma_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE7803_DEG <- SAM_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE7803_DEG <- t_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE7803_DEG$Name)
SAM_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE7803_DEG$Name)
t_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE7803_DEG$Name)

#limma
setwd(datadir)
DEG=read.csv('limma_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% limma_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(limma_GSE7803_DEG)[1]
DEG=merge(DEG,limma_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result=fgsea(pathways = Homo_C2_geneset.list,
                   stats = id,
                   minSize=1,
                   maxSize=10000)
fgsea_result_padj=fgsea_result[order(fgsea_result$padj,decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id,  fgsea_result, gseaParam=0.5)#18*9inches
fgsea_result_nes=fgsea_result[order(abs(fgsea_result$NES),decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id,  fgsea_result, gseaParam=0.5)#18*9inches

setwd(savedir)
limma_nes_top20=fgsea_result_nes$pathway[1:20]
inte=intersect(fgsea_result_padj$pathway[1:20],fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result,file = 'GSE7803_fgsea_limma.rds')
pdf("GSE7803_fgsea_limma.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE7803_fgsea_limma_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[limma_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

#SAM
setwd(datadir)
DEG=read.csv('SAM_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% SAM_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(SAM_GSE7803_DEG)[1]
DEG=merge(DEG,SAM_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$P_adj)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result=fgsea(pathways = Homo_C2_geneset.list,
                   stats = id,
                   minSize=1,
                   maxSize=10000)
fgsea_result_padj=fgsea_result[order(fgsea_result$padj,decreasing = F),]
#plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id,  fgsea_result, gseaParam=0.5)#18*9inches
fgsea_result_nes=fgsea_result[order(abs(fgsea_result$NES),decreasing = T),]
#plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id,  fgsea_result, gseaParam=0.5)#18*9inches

setwd(savedir)
SAM_nes_top20=fgsea_result_nes$pathway[1:20]
inte=intersect(fgsea_result_padj$pathway[1:20],fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result,file = 'GSE7803_fgsea_SAM.rds')
pdf("GSE7803_fgsea_SAM.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE7803_fgsea_SAM_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[SAM_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

#t
setwd(datadir)
DEG=read.csv('t_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% t_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(t_GSE7803_DEG)[1]
DEG=merge(DEG,t_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result=fgsea(pathways = Homo_C2_geneset.list,
                   stats = id,
                   minSize=1,
                   maxSize=10000)
fgsea_result_padj=fgsea_result[order(fgsea_result$padj,decreasing = F),]
#plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id,  fgsea_result, gseaParam=0.5)#18*9inches
fgsea_result_nes=fgsea_result[order(abs(fgsea_result$NES),decreasing = T),]
#plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id,  fgsea_result, gseaParam=0.5)#18*9inches

setwd(savedir)
t_nes_top20=fgsea_result_nes$pathway[1:20]
inte=intersect(fgsea_result_padj$pathway[1:20],fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result,file = 'GSE7803_fgsea_t.rds')
pdf("GSE7803_fgsea_t.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE7803_fgsea_t_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[t_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

commonPathways=intersect(intersect(limma_nes_top20,SAM_nes_top20),t_nes_top20)
saveRDS(commonPathways,file = 'GSE7803_NES_common.rds')


############################################################ GSE9750 fgsea ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
limma_GSE9750_DEG = read.table('limma_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE9750_DEG = read.table('SAM_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE9750_DEG = read.table('t_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE9750_DEG <- limma_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE9750_DEG <- SAM_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE9750_DEG <- t_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE9750_DEG$Name)
SAM_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE9750_DEG$Name)
t_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE9750_DEG$Name)

# limma
setwd(datadir)
DEG=read.csv('limma_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% limma_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(limma_GSE9750_DEG)[1]
DEG=merge(DEG,limma_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result=fgsea(pathways = Homo_C2_geneset.list,
                   stats = id,
                   minSize=1,
                   maxSize=10000)
fgsea_result_padj=fgsea_result[order(fgsea_result$padj,decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id,  fgsea_result, gseaParam=0.5) # 18*9inches
fgsea_result_nes=fgsea_result[order(abs(fgsea_result$NES),decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id,  fgsea_result, gseaParam=0.5) # 18*9inches

setwd(savedir)
limma_nes_top20=fgsea_result_nes$pathway[1:20]
inte=intersect(fgsea_result_padj$pathway[1:20],fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result,file = 'GSE9750_fgsea_limma.rds')
pdf("GSE9750_fgsea_limma.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE9750_fgsea_limma_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[limma_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# SAM
setwd(datadir)
DEG=read.csv('SAM_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% SAM_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(SAM_GSE9750_DEG)[1]
DEG=merge(DEG,SAM_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$P_adj)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result=fgsea(pathways = Homo_C2_geneset.list,
                   stats = id,
                   minSize=1,
                   maxSize=10000)
fgsea_result_padj=fgsea_result[order(fgsea_result$padj,decreasing = F),]
# plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id,  fgsea_result, gseaParam=0.5) # 18*9inches
fgsea_result_nes=fgsea_result[order(abs(fgsea_result$NES),decreasing = T),]
# plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id,  fgsea_result, gseaParam=0.5) # 18*9inches

setwd(savedir)
SAM_nes_top20=fgsea_result_nes$pathway[1:20]
inte=intersect(fgsea_result_padj$pathway[1:20],fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result,file = 'GSE9750_fgsea_SAM.rds')
pdf("GSE9750_fgsea_SAM.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE9750_fgsea_SAM_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[SAM_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# t
setwd(datadir)
DEG=read.csv('t_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% t_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(t_GSE9750_DEG)[1]
DEG=merge(DEG,t_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result=fgsea(pathways = Homo_C2_geneset.list,
                   stats = id,
                   minSize=1,
                   maxSize=10000)
fgsea_result_padj=fgsea_result[order(fgsea_result$padj,decreasing = F),]
# plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id,  fgsea_result, gseaParam=0.5) # 18*9inches
fgsea_result_nes=fgsea_result[order(abs(fgsea_result$NES),decreasing = T),]
# plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id,  fgsea_result, gseaParam=0.5) # 18*9inches

setwd(savedir)
t_nes_top20=fgsea_result_nes$pathway[1:20]
inte=intersect(fgsea_result_padj$pathway[1:20],fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result,file = 'GSE9750_fgsea_t.rds')
pdf("GSE9750_fgsea_t.pdf", width = 17, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE9750_fgsea_t_NES_20.pdf", width = 17, height = 10)
plotGseaTable(Homo_C2_geneset.list[t_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

commonPathways = intersect(intersect(limma_nes_top20, SAM_nes_top20), t_nes_top20)
saveRDS(commonPathways, file = 'GSE9750_NES_common.rds')


############################################################ GSE168244 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

setwd(datadir)
DESeq2_DEG = read.csv('DESeq2_GSE168244_DEG.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE168244_DEG.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE168244_DEG.csv', header = TRUE, row.names = 1)

# DESeq2
DEG = DESeq2_DEG
id = -log(DEG$padj)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
DESeq2_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE168244_fgsea_DESeq2.rds')
pdf("GSE168244_fgsea_DESeq2.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE168244_fgsea_DESeq2_NES_20.pdf", width = 22, height = 10)
plotGseaTable(Homo_C2_geneset.list[DESeq2_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# edgeR
DEG = edgeR_DEG
id = -log(DEG$adj_pval)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
# plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
# plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
edgeR_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE168244_fgsea_edgeR.rds')
pdf("GSE168244_fgsea_edgeR.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE168244_fgsea_edgeR_NES_20.pdf", width = 16, height = 10)
plotGseaTable(Homo_C2_geneset.list[edgeR_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# NOISeq
DEG = NOISeq_DEG
id = DEG$log2FC
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
# plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
# plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
NOISeq_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE168244_fgsea_NOISeq.rds')
pdf("GSE168244_fgsea_NOISeq.pdf", width = 17, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE168244_fgsea_NOISeq_NES_20.pdf", width = 17, height = 10)
plotGseaTable(Homo_C2_geneset.list[NOISeq_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

commonPathways = intersect(intersect(DESeq2_nes_top20, edgeR_nes_top20), NOISeq_nes_top20)
saveRDS(commonPathways, file = 'GSE168244_NES_common.rds')


############################################################ GSE149450 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
DESeq2_DEG = read.csv('DESeq2_GSE149450_DEG_symbol.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE149450_DEG_symbol.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE149450_DEG_symbol.csv', header = TRUE, row.names = 1)

DESeq2_DEG=DESeq2_DEG %>% filter(symbol != '')
edgeR_DEG=edgeR_DEG %>% filter(symbol != '')
NOISeq_DEG=NOISeq_DEG %>% filter(symbol != '')

DESeq2_DEG$X=rownames(DESeq2_DEG)
edgeR_DEG$X=rownames(edgeR_DEG)
NOISeq_DEG$X=rownames(NOISeq_DEG)

# DESeq2
setwd(datadir)
DEG = read.csv('DESeq2_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$padj)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
DESeq2_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE149450_fgsea_DESeq2.rds')
pdf("GSE149450_fgsea_DESeq2.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE149450_fgsea_DESeq2_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[DESeq2_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# edgeR
setwd(datadir)
DEG = read.csv('edgeR_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$adj_pval)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
edgeR_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE149450_fgsea_edgeR.rds')
pdf("GSE149450_fgsea_edgeR.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE149450_fgsea_edgeR_NES_20.pdf", width = 16, height = 10)
plotGseaTable(Homo_C2_geneset.list[edgeR_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# NOISeq
setwd(datadir)
DEG = read.csv('NOISeq_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% NOISeq_DEG$X,]
DEG = merge(DEG, NOISeq_DEG, by = colnames(DEG)[1])
id = DEG$log2FC
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
# plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
# plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
NOISeq_nes_top20 = fgsea_result_nes$pathway[1:20]
NOISeq_nes_top30 = fgsea_result_nes$pathway[1:30]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE149450_fgsea_NOISeq.rds')
pdf("GSE149450_fgsea_NOISeq.pdf", width = 16, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE149450_fgsea_NOISeq_NES_20.pdf", width = 16, height = 10)
plotGseaTable(Homo_C2_geneset.list[NOISeq_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

commonPathways = intersect(intersect(DESeq2_nes_top20, edgeR_nes_top20), NOISeq_nes_top20)
saveRDS(commonPathways, file = 'GSE149450_NES_common.rds')

#### enrichment of the common in GSE7803 PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
limma_GSE7803_DEG = read.table('limma_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE7803_DEG = read.table('SAM_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE7803_DEG = read.table('t_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE7803_DEG <- limma_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE7803_DEG <- SAM_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE7803_DEG <- t_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE7803_DEG$Name)
SAM_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE7803_DEG$Name)
t_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE7803_DEG$Name)

#limma
setwd(datadir)
DEG=read.csv('limma_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% limma_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(limma_GSE7803_DEG)[1]
DEG=merge(DEG,limma_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
setwd(savedir)
pdf('GSE7803_limma_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +                id) + labs(title="GSE7803_limma_fgsea")
dev.off()

#SAM
setwd(datadir)
DEG=read.csv('SAM_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% SAM_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(SAM_GSE7803_DEG)[1]
DEG=merge(DEG,SAM_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$P_adj)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE7803_SAM_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +                id) + labs(title="GSE7803_SAM_fgsea")
dev.off()


#t
setwd(datadir)
DEG=read.csv('t_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% t_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(t_GSE7803_DEG)[1]
DEG=merge(DEG,t_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE7803_t_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +                id) + labs(title="GSE7803_t_fgsea")
dev.off()


########################################## enrichment of the common in GSE9750 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
limma_GSE9750_DEG = read.table('limma_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE9750_DEG = read.table('SAM_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE9750_DEG = read.table('t_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE9750_DEG <- limma_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE9750_DEG <- SAM_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE9750_DEG <- t_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE9750_DEG$Name)
SAM_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE9750_DEG$Name)
t_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE9750_DEG$Name)

# limma
setwd(datadir)
DEG=read.csv('limma_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% limma_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(limma_GSE9750_DEG)[1]
DEG=merge(DEG,limma_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE9750_limma_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +                id) + labs(title="GSE9750_limma_fgsea")
dev.off()

# SAM
setwd(datadir)
DEG=read.csv('SAM_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% SAM_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(SAM_GSE9750_DEG)[1]
DEG=merge(DEG,SAM_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$P_adj)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE9750_SAM_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +                id) + labs(title="GSE9750_SAM_fgsea")
dev.off()

# t
setwd(datadir)
DEG=read.csv('t_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% t_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(t_GSE9750_DEG)[1]
DEG=merge(DEG,t_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE9750_t_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +                id) + labs(title="GSE9750_t_fgsea")
dev.off()

########################################## enrichment of the common in GSE168244 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

setwd(datadir)
DESeq2_DEG = read.csv('DESeq2_GSE168244_DEG.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE168244_DEG.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE168244_DEG.csv', header = TRUE, row.names = 1)

# DESeq2
DEG = DESeq2_DEG
id = -log(DEG$padj)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE168244_DESeq2_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
                 +id) + labs(title="GSE168244_DESeq2_fgsea")+
                               theme(plot.title = element_text(hjust = 1))
dev.off()

# edgeR
DEG = edgeR_DEG
id = -log(DEG$adj_pval)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE168244_edgeR_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +id) + labs(title="GSE168244_edgeR_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

# NOISeq
DEG = NOISeq_DEG
id = DEG$log2FC
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE168244_NOISeq_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               + id) + labs(title="GSE168244_NOISeq_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

########################################## enrichment of the common in GSE149450 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
DESeq2_DEG = read.csv('DESeq2_GSE149450_DEG_symbol.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE149450_DEG_symbol.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE149450_DEG_symbol.csv', header = TRUE, row.names = 1)

DESeq2_DEG=DESeq2_DEG %>% filter(symbol != '')
edgeR_DEG=edgeR_DEG %>% filter(symbol != '')
NOISeq_DEG=NOISeq_DEG %>% filter(symbol != '')

DESeq2_DEG$X=rownames(DESeq2_DEG)
edgeR_DEG$X=rownames(edgeR_DEG)
NOISeq_DEG$X=rownames(NOISeq_DEG)

# DESeq2
setwd(datadir)
DEG = read.csv('DESeq2_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$padj)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE149450_DESeq2_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +               id) + labs(title="GSE149450_DESeq2_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

# edgeR
setwd(datadir)
DEG = read.csv('edgeR_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$adj_pval)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE149450_edgeR_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
               +               id) + labs(title="GSE149450_edgeR_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()


# NOISeq
setwd(datadir)
DEG = read.csv('NOISeq_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% NOISeq_DEG$X,]
DEG = merge(DEG, NOISeq_DEG, by = colnames(DEG)[1])
id = DEG$log2FC
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE149450_NOISeq_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP"]],
                 +               id) + labs(title="GSE149450_NOISeq_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

#### enrichment of the common in GSE7803 ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
limma_GSE7803_DEG = read.table('limma_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE7803_DEG = read.table('SAM_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE7803_DEG = read.table('t_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE7803_DEG <- limma_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE7803_DEG <- SAM_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE7803_DEG <- t_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE7803_DEG$Name)
SAM_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE7803_DEG$Name)
t_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE7803_DEG$Name)

#limma
setwd(datadir)
DEG=read.csv('limma_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% limma_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(limma_GSE7803_DEG)[1]
DEG=merge(DEG,limma_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

library(fgsea)
library(ggplot2)
setwd(savedir)
pdf('GSE7803_limma_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +                id) + labs(title="GSE7803_limma_fgsea")
dev.off()

#SAM
setwd(datadir)
DEG=read.csv('SAM_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% SAM_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(SAM_GSE7803_DEG)[1]
DEG=merge(DEG,SAM_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$P_adj)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE7803_SAM_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +                id) + labs(title="GSE7803_SAM_fgsea")
dev.off()


#t
setwd(datadir)
DEG=read.csv('t_GSE7803_DEG.csv')
DEG = DEG[DEG$X %in% t_GSE7803_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(t_GSE7803_DEG)[1]
DEG=merge(DEG,t_GSE7803_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE7803_t_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +                id) + labs(title="GSE7803_t_fgsea")
dev.off()
########################################## enrichment of the common in GSE9750 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
limma_GSE9750_DEG = read.table('limma_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE9750_DEG = read.table('SAM_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE9750_DEG = read.table('t_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE9750_DEG <- limma_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE9750_DEG <- SAM_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE9750_DEG <- t_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE9750_DEG$Name)
SAM_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE9750_DEG$Name)
t_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE9750_DEG$Name)

# limma
setwd(datadir)
DEG=read.csv('limma_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% limma_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(limma_GSE9750_DEG)[1]
DEG=merge(DEG,limma_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE9750_limma_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +                id) + labs(title="GSE9750_limma_fgsea")
dev.off()

# SAM
setwd(datadir)
DEG=read.csv('SAM_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% SAM_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(SAM_GSE9750_DEG)[1]
DEG=merge(DEG,SAM_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$P_adj)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE9750_SAM_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +                id) + labs(title="GSE9750_SAM_fgsea")
dev.off()

# t
setwd(datadir)
DEG=read.csv('t_GSE9750_DEG.csv')
DEG = DEG[DEG$X %in% t_GSE9750_DEG$AFFYMETRIX_3PRIME_IVT_ID, ]
colnames(DEG)[1]=colnames(t_GSE9750_DEG)[1]
DEG=merge(DEG,t_GSE9750_DEG,by=colnames(DEG)[1])
id=-log(DEG$adj_pval)
names(id)=DEG$symbol
id=id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE9750_t_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +                id) + labs(title="GSE9750_t_fgsea")
dev.off()

########################################## enrichment of the common in GSE168244 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

setwd(datadir)
DESeq2_DEG = read.csv('DESeq2_GSE168244_DEG.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE168244_DEG.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE168244_DEG.csv', header = TRUE, row.names = 1)

# DESeq2
DEG = DESeq2_DEG
id = -log(DEG$padj)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE168244_DESeq2_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +id) + labs(title="GSE168244_DESeq2_fgsea")+
  theme(plot.title = element_text(hjust = 1))
dev.off()

# edgeR
DEG = edgeR_DEG
id = -log(DEG$adj_pval)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE168244_edgeR_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +id) + labs(title="GSE168244_edgeR_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

# NOISeq
DEG = NOISeq_DEG
id = DEG$log2FC
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE168244_NOISeq_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               + id) + labs(title="GSE168244_NOISeq_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

########################################## enrichment of the common in GSE149450 ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
DESeq2_DEG = read.csv('DESeq2_GSE149450_DEG_symbol.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE149450_DEG_symbol.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE149450_DEG_symbol.csv', header = TRUE, row.names = 1)

DESeq2_DEG=DESeq2_DEG %>% filter(symbol != '')
edgeR_DEG=edgeR_DEG %>% filter(symbol != '')
NOISeq_DEG=NOISeq_DEG %>% filter(symbol != '')

DESeq2_DEG$X=rownames(DESeq2_DEG)
edgeR_DEG$X=rownames(edgeR_DEG)
NOISeq_DEG$X=rownames(NOISeq_DEG)

# DESeq2
setwd(datadir)
DEG = read.csv('DESeq2_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$padj)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE149450_DESeq2_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +               id) + labs(title="GSE149450_DESeq2_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()

# edgeR
setwd(datadir)
DEG = read.csv('edgeR_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$adj_pval)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE149450_edgeR_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +               id) + labs(title="GSE149450_edgeR_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()


# NOISeq
setwd(datadir)
DEG = read.csv('NOISeq_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% NOISeq_DEG$X,]
DEG = merge(DEG, NOISeq_DEG, by = colnames(DEG)[1])
id = DEG$log2FC
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

setwd(savedir)
pdf('GSE149450_NOISeq_fgsea.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               +               id) + labs(title="GSE149450_NOISeq_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()


#### shared 22 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE149450")
GSE149450_fgsea_NOISeq=readRDS('GSE149450_fgsea_NOISeq.rds')
GSE149450_fgsea_edgeR=readRDS('GSE149450_fgsea_edgeR.rds')
GSE149450_fgsea_DESeq2=readRDS('GSE149450_fgsea_DESeq2.rds')
vector=c("GSE149450_fgsea_NOISeq","GSE149450_fgsea_edgeR","GSE149450_fgsea_DESeq2")
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$pval <0.05),]
  res=res[,c(1,6)]
  colnames(res)[2]=vector[i]
  res=res[order(abs(res[,2])),]
  if(i==1){
    res_all=res
  }else{
    res_all=merge(res_all,res,by="pathway")
  }
}

res_all$mean_nes=rowMeans(res_all[,2:4])
GSE149450_enrich=res_all[which(abs(res_all$mean_nes)>1),]

saveRDS(GSE149450_enrich,file = 'GSE149450_enrich.rds')
#### shared 203  ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE168244")
GSE168244_fgsea_NOISeq=readRDS('GSE168244_fgsea_NOISeq.rds')
GSE168244_fgsea_edgeR=readRDS('GSE168244_fgsea_edgeR.rds')
GSE168244_fgsea_DESeq2=readRDS('GSE168244_fgsea_DESeq2.rds')
vector=c("GSE168244_fgsea_NOISeq","GSE168244_fgsea_edgeR","GSE168244_fgsea_DESeq2")
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$pval <0.05),]
  res=res[,c(1,6)]
  colnames(res)[2]=vector[i]
  res=res[order(abs(res[,2])),]
  if(i==1){
    res_all=res
  }else{
    res_all=merge(res_all,res,by="pathway")
  }
}

res_all$mean_nes=rowMeans(res_all[,2:4])
GSE168244_enrich=res_all[which(abs(res_all$mean_nes)>1),]
saveRDS(GSE168244_enrich,file = 'GSE168244_enrich.rds')
#### shared 735 ####
rm(list = ls())
gc()
wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE9750"
setwd(wd)
dir(pattern = '.rds',path = wd)
GSE9750_fgsea_limma=readRDS("GSE9750_fgsea_limma.rds")
GSE9750_fgsea_SAM=readRDS("GSE9750_fgsea_SAM.rds")
GSE9750_fgsea_t=readRDS("GSE9750_fgsea_t.rds" )
vector=c("GSE9750_fgsea_limma","GSE9750_fgsea_SAM","GSE9750_fgsea_t" )
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$pval <0.05),]
  res=res[,c(1,6)]
  colnames(res)[2]=vector[i]
  res=res[order(abs(res[,2])),]
  if(i==1){
    res_all=res
  }else{
    res_all=merge(res_all,res,by="pathway")
  }
}

res_all$mean_nes=rowMeans(res_all[,2:4])
GSE9750_enrich=res_all[which(abs(res_all$mean_nes)>1),]
saveRDS(GSE9750_enrich,file = 'GSE9750_enrich.rds')

#### shared 596 ####
rm(list = ls())
gc()
wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE7803"
setwd(wd)
dir(pattern = '.rds',path = wd)
GSE7803_fgsea_limma=readRDS('GSE7803_fgsea_limma.rds')
GSE7803_fgsea_SAM=readRDS('GSE7803_fgsea_SAM.rds')
GSE7803_fgsea_t=readRDS('GSE7803_fgsea_t.rds')
vector=c('GSE7803_fgsea_limma','GSE7803_fgsea_SAM','GSE7803_fgsea_t')
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$pval <0.05),]
  res=res[,c(1,6)]
  colnames(res)[2]=vector[i]
  res=res[order(abs(res[,2])),]
  if(i==1){
    res_all=res
  }else{
    res_all=merge(res_all,res,by="pathway")
  }
}

res_all$mean_nes=rowMeans(res_all[,2:4])
GSE7803_enrich=res_all[which(abs(res_all$mean_nes)>1),]
saveRDS(GSE7803_enrich,file = 'GSE7803_enrich.rds')

######################## save and collect results -- common & C 4 3 ####



testx=intersect(GSE149450_enrich$pathway,GSE168244_enrich$pathway)
testx1=intersect(GSE149450_enrich$pathway,GSE7803_enrich$pathway)
testx2=intersect(GSE149450_enrich$pathway,GSE9750_enrich$pathway)
testy=intersect(GSE7803_enrich$pathway,GSE9750_enrich$pathway)
common=intersect(testx,testy)

GSE149450_enrich$abs=abs(GSE149450_enrich$mean_nes)# 12th / 22
GSE168244_enrich$abs=abs(GSE168244_enrich$mean_nes)# 34th / 203
GSE9750_enrich$abs=abs(GSE9750_enrich$mean_nes)# 444th / 735
GSE7803_enrich$abs=abs(GSE7803_enrich$mean_nes)# 132th /596


df=data.frame()
df=GSE149450_fgsea_DESeq2[3006,]
df=rbind(df,GSE149450_fgsea_edgeR[2994,])
df=rbind(df,GSE149450_fgsea_NOISeq[3032,])
df=rbind(df,GSE168244_fgsea_DESeq2[2999,])
df=rbind(df,GSE168244_fgsea_edgeR[2996,])
df=rbind(df,GSE168244_fgsea_NOISeq[3037,])
df=rbind(df,GSE7803_fgsea_limma[2991,])
df=rbind(df,GSE7803_fgsea_SAM[2988,])
df=rbind(df,GSE7803_fgsea_t[2989,])
df=rbind(df,GSE9750_fgsea_limma[3043,])
df=rbind(df,GSE9750_fgsea_SAM[3036,])
df=rbind(df,GSE9750_fgsea_t[3035,])

df$inf=c('GSE149450_fgsea_DESeq2','GSE149450_fgsea_edgeR','GSE149450_fgsea_NOISeq'
         ,'GSE168244_fgsea_DESeq2','GSE168244_fgsea_edgeR','GSE168244_fgsea_NOISeq'
         ,'GSE7803_fgsea_limma','GSE7803_fgsea_SAM','GSE7803_fgsea_t'
         ,'GSE9750_fgsea_limma','GSE9750_fgsea_SAM','GSE9750_fgsea_t')

wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
setwd(wd)
saveRDS(df,file = 'result_df.rds')
write.csv(df[,-8],file = 'result_df.csv',row.names=F)


df_C_4_3=rbind(GSE149450_enrich[,c(1,5,6)],GSE168244_enrich[,c(1,5,6)])
df_C_4_3=rbind(df_C_4_3,GSE9750_enrich[,c(1,5,6)])
df_C_4_3=rbind(df_C_4_3,GSE7803_enrich[,c(1,5,6)])
df_C_4_3_gene=data.frame()
df_C_4_3_gene=table(df_C_4_3$pathway)
result=subset(df_C_4_3_gene, df_C_4_3_gene > 2)
result_table=as.data.frame(result)
saveRDS(result_table,file = 'result_C_4_3.rds')


vec=c('GSE149450_enrich','GSE168244_enrich','GSE9750_enrich','GSE7803_enrich')
for (i in 1:length(vec)){
  res=get(vec[i])
  res=res[res$pathway%in% result_table$Var1,]
  res$State=ifelse(res$mean_nes>0,'up','down')
  res=res[,c(1,6,7)]
  t=nrow(res)
  res=rbind(res,res,res)
  if(i<3){res$Method=c(rep("NOISeq", t), rep("edgeR", t), rep("DESeq2", t))
  }else{res$Method=c(rep("limma", t), rep("SAM", t), rep("t", t))}
  sample_name=sub("_enrich$", "", vec[i])
  res$Sample=rep(sample_name,nrow(res))
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
}
colnames(res_all)[1]='Pathway'
final=res_all[,c('Sample','Method','Pathway','State')]
saveRDS(res_all,file = 'Sankey_data_C_4_3_NES.rds')
write.csv(final,file = 'Sankey_data_C_4_3.csv',row.names = F)











#### GSE149450 NES >3 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE149450")
GSE149450_fgsea_NOISeq=readRDS('GSE149450_fgsea_NOISeq.rds')
GSE149450_fgsea_edgeR=readRDS('GSE149450_fgsea_edgeR.rds')
GSE149450_fgsea_DESeq2=readRDS('GSE149450_fgsea_DESeq2.rds')
vector=c("GSE149450_fgsea_NOISeq","GSE149450_fgsea_edgeR","GSE149450_fgsea_DESeq2")
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>3),]
  if(i==1){res$Method=rep('NOISeq',nrow(res))
  }else if(i==2){res$Method=rep('edgeR',nrow(res))
  }else{res$Method=rep('DESeq2',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
saveRDS(res_all,file = 'GSE149450_NES_3.rds')
#### GSE168244 NES >3 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE168244")
GSE168244_fgsea_NOISeq=readRDS('GSE168244_fgsea_NOISeq.rds')
GSE168244_fgsea_edgeR=readRDS('GSE168244_fgsea_edgeR.rds')
GSE168244_fgsea_DESeq2=readRDS('GSE168244_fgsea_DESeq2.rds')
vector=c("GSE168244_fgsea_NOISeq","GSE168244_fgsea_edgeR","GSE168244_fgsea_DESeq2")
res_all=data.frame()

for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>3),]
  if(i==1){res$Method=rep('NOISeq',nrow(res))
  }else if(i==2){res$Method=rep('edgeR',nrow(res))
  }else{res$Method=rep('DESeq2',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
saveRDS(res_all,file = 'GSE168244_NES_3.rds')

#### GSE9750 NES >3 #######
wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE9750"
setwd(wd)
dir(pattern = '.rds',path = wd)
GSE9750_fgsea_limma=readRDS("GSE9750_fgsea_limma.rds")
GSE9750_fgsea_SAM=readRDS("GSE9750_fgsea_SAM.rds")
GSE9750_fgsea_t=readRDS("GSE9750_fgsea_t.rds" )
vector=c("GSE9750_fgsea_limma","GSE9750_fgsea_SAM","GSE9750_fgsea_t" )
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>3),]
  if(i==1){res$Method=rep('limma',nrow(res))
  }else if(i==2){res$Method=rep('SAM',nrow(res))
  }else{res$Method=rep('t',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
res_all$Sample=rep('GSE9750',nrow(res_all))
saveRDS(res_all,file = 'GSE9750_NES_3.rds')
#### GSE7803 NES >3 ###########
wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE7803"
setwd(wd)
dir(pattern = '.rds',path = wd)
GSE7803_fgsea_limma=readRDS('GSE7803_fgsea_limma.rds')
GSE7803_fgsea_SAM=readRDS('GSE7803_fgsea_SAM.rds')
GSE7803_fgsea_t=readRDS('GSE7803_fgsea_t.rds')
vector=c('GSE7803_fgsea_limma','GSE7803_fgsea_SAM','GSE7803_fgsea_t')
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>3),]
  if(i==1){res$Method=rep('limma',nrow(res))
  }else if(i==2){res$Method=rep('SAM',nrow(res))
  }else{res$Method=rep('t',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
res_all$Sample=rep('GSE7803',nrow(res_all))
saveRDS(res_all,file = 'GSE7803_NES_3.rds')
######################## save and collect results -- NES >3 & NES >3.09 ####
gc()
rm(list = ls())
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE149450")
GSE149450=readRDS('GSE149450_NES_3.rds')
GSE149450$Sample=rep('GSE149450',nrow(GSE149450))
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE168244")
GSE168244=readRDS('GSE168244_NES_3.rds')
GSE168244$Sample=rep('GSE168244',nrow(GSE168244))
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE9750")
GSE9750=readRDS('GSE9750_NES_3.rds')
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE7803")
GSE7803=readRDS('GSE7803_NES_3.rds')
result_NES_3=rbind(GSE149450,GSE168244,GSE9750,GSE7803)
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
saveRDS(result_NES_3,file = 'result_NES_3.rds')
colnames(result_NES_3)[1]='Pathway'
write.csv(result_NES_3[,c('Sample','Method','Pathway')],file = 'result_NES_3.csv',row.names = F)


result_NES_309=result_NES_3[which(abs(result_NES_3$NES)>3.09),]
write.csv(result_NES_309[,c('Sample','Method','Pathway')],file = 'result_NES_309.csv',row.names = F)

######################## add connection with key genes #### 
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
result_NES_309 <- readRDS("~/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/result_NES_3.rds")
KG=c('CDK1','TOP2A','BRCA1','CDC20','FEN1','BUB1B','CDC45','KIF23','KIF11','MCM7')
result_NES_309=result_NES_309[which(abs(result_NES_3$NES)>3.09),]
Sankey=result_NES_309[,c(1,9,10)]

filterpathway <- c()
for (i in 1:nrow(Sankey)){
  show=intersect(KG,unlist(result_NES_309[i, 8]))
  if(length(show)>0){
    add=Sankey[i,1]
    if(i==1){
      filterpathway=add
    }else{
      filterpathway=c(filterpathway,add)
    }}
}

Sankey=result_NES_309[which(result_NES_309$pathway%in% filterpathway),]
Sankey_data=data.frame()
for (i in 1:nrow(Sankey)){
  show=intersect(KG,unlist(Sankey[i, 8]))
  new_rows=do.call(rbind,replicate(length(show),Sankey[i,],simplify = F))
  new_rows$KeyGene=show
  if(i==1){
    Sankey_data=new_rows
  }else{
    Sankey_data=rbind(Sankey_data,new_rows)
  }
}
colnames(Sankey_data)[1]='Pathway'
saveRDS(Sankey_data,file = 'Sankey_data_NES_3_09_KG.rds')
write.csv(Sankey_data[,c('Sample','Method','Pathway','KeyGene')],file = 'Sankey_data_NES_3_09_KG.csv',row.names = F)





#### GSE149450 NES >2 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE149450")
GSE149450_fgsea_NOISeq=readRDS('GSE149450_fgsea_NOISeq.rds')
GSE149450_fgsea_edgeR=readRDS('GSE149450_fgsea_edgeR.rds')
GSE149450_fgsea_DESeq2=readRDS('GSE149450_fgsea_DESeq2.rds')
vector=c("GSE149450_fgsea_NOISeq","GSE149450_fgsea_edgeR","GSE149450_fgsea_DESeq2")
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>2),]
  if(i==1){res$Method=rep('NOISeq',nrow(res))
  }else if(i==2){res$Method=rep('edgeR',nrow(res))
  }else{res$Method=rep('DESeq2',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
saveRDS(res_all,file = 'GSE149450_NES_2.rds')
#### GSE168244 NES >2 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE168244")
GSE168244_fgsea_NOISeq=readRDS('GSE168244_fgsea_NOISeq.rds')
GSE168244_fgsea_edgeR=readRDS('GSE168244_fgsea_edgeR.rds')
GSE168244_fgsea_DESeq2=readRDS('GSE168244_fgsea_DESeq2.rds')
vector=c("GSE168244_fgsea_NOISeq","GSE168244_fgsea_edgeR","GSE168244_fgsea_DESeq2")
res_all=data.frame()

for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>2),]
  if(i==1){res$Method=rep('NOISeq',nrow(res))
  }else if(i==2){res$Method=rep('edgeR',nrow(res))
  }else{res$Method=rep('DESeq2',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
saveRDS(res_all,file = 'GSE168244_NES_2.rds')

#### GSE9750 NES >2 #######
wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE9750"
setwd(wd)
dir(pattern = '.rds',path = wd)
GSE9750_fgsea_limma=readRDS("GSE9750_fgsea_limma.rds")
GSE9750_fgsea_SAM=readRDS("GSE9750_fgsea_SAM.rds")
GSE9750_fgsea_t=readRDS("GSE9750_fgsea_t.rds" )
vector=c("GSE9750_fgsea_limma","GSE9750_fgsea_SAM","GSE9750_fgsea_t" )
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>2),]
  if(i==1){res$Method=rep('limma',nrow(res))
  }else if(i==2){res$Method=rep('SAM',nrow(res))
  }else{res$Method=rep('t',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
res_all$Sample=rep('GSE9750',nrow(res_all))
saveRDS(res_all,file = 'GSE9750_NES_2.rds')
#### GSE7803 NES >2 ###########
wd="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE7803"
setwd(wd)
dir(pattern = '.rds',path = wd)
GSE7803_fgsea_limma=readRDS('GSE7803_fgsea_limma.rds')
GSE7803_fgsea_SAM=readRDS('GSE7803_fgsea_SAM.rds')
GSE7803_fgsea_t=readRDS('GSE7803_fgsea_t.rds')
vector=c('GSE7803_fgsea_limma','GSE7803_fgsea_SAM','GSE7803_fgsea_t')
res_all=data.frame()
for (i in 1:length(vector)){
  res=get(vector[i])
  res=res[which(res$padj <0.01),]
  res=res[which(abs(res$NES)>2),]
  if(i==1){res$Method=rep('limma',nrow(res))
  }else if(i==2){res$Method=rep('SAM',nrow(res))
  }else{res$Method=rep('t',nrow(res))}
  if(i==1){
    res_all=res
  }else{
    res_all=rbind(res_all,res)
  }
  
}
res_all$Sample=rep('GSE7803',nrow(res_all))
saveRDS(res_all,file = 'GSE7803_NES_2.rds')
######################## save and collect results -- NES >2  ####
gc()
rm(list = ls())
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE149450")
GSE149450=readRDS('GSE149450_NES_2.rds')
GSE149450$Sample=rep('GSE149450',nrow(GSE149450))
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE168244")
GSE168244=readRDS('GSE168244_NES_2.rds')
GSE168244$Sample=rep('GSE168244',nrow(GSE168244))
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE9750")
GSE9750=readRDS('GSE9750_NES_2.rds')
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE7803")
GSE7803=readRDS('GSE7803_NES_2.rds')
result_NES_2=rbind(GSE149450,GSE168244,GSE9750,GSE7803)
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
colnames(result_NES_2)[1]='Pathway'
saveRDS(result_NES_2,file = 'result_NES_2.rds')
write.csv(result_NES_2[,c('Sample','Method','Pathway')],file = 'result_NES_2.csv',row.names = F)


######################## add connection with key genes #### 
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
result_NES_2 <- readRDS("~/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/result_NES_2.rds")
KG=c('CDK1','TOP2A','BRCA1','CDC20','FEN1','BUB1B','CDC45','KIF23','KIF11','MCM7')
Sankey=result_NES_2[,c(1,9,10)]

filterpathway <- c()
for (i in 1:nrow(Sankey)){
  show=intersect(KG,unlist(result_NES_2[i, 8]))
  if(length(show)>0){
    add=Sankey[i,1]
    if(i==1){
      filterpathway=add
    }else{
      filterpathway=c(filterpathway,add)
    }}
}

Sankey=result_NES_2[which(result_NES_2$Pathway%in% filterpathway),]
Sankey_data=data.frame()
for (i in 1:nrow(Sankey)){
  show=intersect(KG,unlist(Sankey[i, 8]))
  new_rows=do.call(rbind,replicate(length(show),Sankey[i,],simplify = F))
  new_rows$KeyGene=show
  if(i==1){
    Sankey_data=new_rows
  }else{
    Sankey_data=rbind(Sankey_data,new_rows)
  }
}
saveRDS(Sankey_data,file = 'Sankey_data_NES_2_KG.rds')
write.csv(Sankey_data[,c('Sample','Method','Pathway','KeyGene')],file = 'Sankey_data_NES_2_KG.csv',row.names = F)



############################################################ GSE168244 fgsea theta ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

setwd(datadir)
DESeq2_DEG = read.csv('DESeq2_GSE168244_DEG.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE168244_DEG.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE168244_DEG.csv', header = TRUE, row.names = 1)

# DESeq2
DEG = DESeq2_DEG
id = -log(DEG$padj)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
DESeq2_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE168244_fgsea_DESeq2.rds')
pdf("GSE168244_fgsea_DESeq2.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE168244_fgsea_DESeq2_NES_20.pdf", width = 22, height = 10)
plotGseaTable(Homo_C2_geneset.list[DESeq2_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# edgeR
DEG = edgeR_DEG
id = -log(DEG$adj_pval)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
# plotGseaTable(Homo_all_geneset[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
# plotGseaTable(Homo_all_geneset[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
edgeR_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE168244_fgsea_edgeR.rds')
pdf("GSE168244_fgsea_edgeR.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE168244_fgsea_edgeR_NES_20.pdf", width = 16, height = 10)
plotGseaTable(Homo_C2_geneset.list[edgeR_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# NOISeq_abs_theta
DEG = NOISeq_DEG
id = abs(DEG$theta)
names(id) = rownames(DEG)
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
saveRDS(fgsea_result, file = 'GSE168244_fgsea_NOISeq_theta_abs.rds')

pdf('GSE168244_NOISeq_fgsea_theta_abs.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               + id) + labs(title="GSE168244_NOISeq_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()



############################################################ GSE149450 fgsea theta ####
rm(list = ls())
gc()
datadir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/DE'
savedir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea"
IDdir="/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/IDconversion"

library(msigdbr)
Homo_C2_geneset = msigdbr(species = "Homo sapiens",
                          category = "C2") %>% dplyr::select(gs_name,gene_symbol)
Homo_C2_geneset.list <- Homo_C2_geneset %>% split(x = .$gene_symbol, f = .$gs_name)

library(dplyr)
setwd(IDdir)
DESeq2_DEG = read.csv('DESeq2_GSE149450_DEG_symbol.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE149450_DEG_symbol.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE149450_DEG_symbol.csv', header = TRUE, row.names = 1)

DESeq2_DEG=DESeq2_DEG %>% filter(symbol != '')
edgeR_DEG=edgeR_DEG %>% filter(symbol != '')
NOISeq_DEG=NOISeq_DEG %>% filter(symbol != '')

DESeq2_DEG$X=rownames(DESeq2_DEG)
edgeR_DEG$X=rownames(edgeR_DEG)
NOISeq_DEG$X=rownames(NOISeq_DEG)

# DESeq2
setwd(datadir)
DEG = read.csv('DESeq2_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$padj)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
DESeq2_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE149450_fgsea_DESeq2.rds')
pdf("GSE149450_fgsea_DESeq2.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE149450_fgsea_DESeq2_NES_20.pdf", width = 12, height = 10)
plotGseaTable(Homo_C2_geneset.list[DESeq2_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# edgeR
setwd(datadir)
DEG = read.csv('edgeR_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% DESeq2_DEG$X,]
DEG = merge(DEG, DESeq2_DEG, by = colnames(DEG)[1])
id = -log(DEG$adj_pval)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)
fgsea_result_padj = fgsea_result[order(fgsea_result$padj, decreasing = F),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_padj$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches
fgsea_result_nes = fgsea_result[order(abs(fgsea_result$NES), decreasing = T),]
plotGseaTable(Homo_C2_geneset.list[fgsea_result_nes$pathway[1:10]], id, fgsea_result, gseaParam = 0.5) # 18*9inches

setwd(savedir)
edgeR_nes_top20 = fgsea_result_nes$pathway[1:20]
inte = intersect(fgsea_result_padj$pathway[1:20], fgsea_result_nes$pathway[1:20])
saveRDS(fgsea_result, file = 'GSE149450_fgsea_edgeR.rds')
pdf("GSE149450_fgsea_edgeR.pdf", width = 12, height = 3)
plotGseaTable(Homo_C2_geneset.list[inte], id, fgsea_result, gseaParam = 0.5)
dev.off()

pdf("GSE149450_fgsea_edgeR_NES_20.pdf", width = 16, height = 10)
plotGseaTable(Homo_C2_geneset.list[edgeR_nes_top20], id, fgsea_result, gseaParam = 0.5)
dev.off()

# NOISeq_theta
setwd(datadir)
DEG = read.csv('NOISeq_GSE149450_DEG.csv')
DEG = DEG[DEG$X %in% NOISeq_DEG$X,]
DEG = merge(DEG, NOISeq_DEG, by = colnames(DEG)[1])
id = abs(DEG$theta)
names(id) = DEG$symbol
id = id[order(id, decreasing = TRUE)]

library(fgsea)
fgsea_result = fgsea(pathways = Homo_C2_geneset.list,
                     stats = id,
                     minSize = 1,
                     maxSize = 10000)

setwd(savedir)
saveRDS(fgsea_result, file = 'GSE149450_fgsea_NOISeq_theta_abs.rds')

pdf('GSE149450_NOISeq_fgsea_theta_abs.pdf',width =2.5,height = 2)
plotEnrichment(Homo_C2_geneset.list[["ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER"]],
               + id) + labs(title="GSE149450_NOISeq_fgsea")+
  theme(plot.title = element_text(hjust = 1)) 
dev.off()






