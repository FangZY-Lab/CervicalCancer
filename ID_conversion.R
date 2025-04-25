#Microarray:Gene ID conversion (https://david.ncifcrf.gov/conversion.jsp)

# GSE7803 #################################
#first use david online

rm(list = ls())
gc()
library(dplyr)
limma_GSE7803_DEG = read.table('limma_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE7803_DEG = read.table('SAM_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE7803_DEG = read.table('t_GSE7803_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE7803_DEG <- limma_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE7803_DEG <- SAM_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE7803_DEG <- t_GSE7803_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE7803_DEG$Name)
SAM_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE7803_DEG$Name)
t_GSE7803_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE7803_DEG$Name)

limma_GSE7803_DEG_rm=unique(limma_GSE7803_DEG$symbol)
SAM_GSE7803_DEG_rm=unique(SAM_GSE7803_DEG$symbol)
t_GSE7803_DEG_rm=unique(t_GSE7803_DEG$symbol)

GSE7803_DEG_shared=intersect(limma_GSE7803_DEG_rm,SAM_GSE7803_DEG_rm)
GSE7803_DEG_shared=intersect(GSE7803_DEG_shared,t_GSE7803_DEG_rm)

write.csv(limma_GSE7803_DEG,file = 'limma_GSE7803_DEG_symbol.csv')
write.csv(SAM_GSE7803_DEG,file = 'SAM_GSE7803_DEG_symbol.csv')
write.csv(t_GSE7803_DEG,file = 't_GSE7803_DEG_symbol.csv')
saveRDS(GSE7803_DEG_shared,file = 'GSE7803_DEG_shared.rds')
write.csv(GSE7803_DEG_shared, file = "GSE7803_DEG_shared.csv", row.names = FALSE)

# GSE9750 #################################
#first use david online

library(dplyr)
limma_GSE9750_DEG = read.table('limma_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
SAM_GSE9750_DEG = read.table('SAM_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")
t_GSE9750_DEG = read.table('t_GSE9750_DEG.txt', sep = "\t", header = TRUE, quote = "")

limma_GSE9750_DEG <- limma_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
SAM_GSE9750_DEG <- SAM_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")
t_GSE9750_DEG <- t_GSE9750_DEG %>% filter(AFFYMETRIX_3PRIME_IVT_ID != "")

limma_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", limma_GSE9750_DEG$Name)
SAM_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", SAM_GSE9750_DEG$Name)
t_GSE9750_DEG$symbol=sub(".*\\(([^)]+)\\).*", "\\1", t_GSE9750_DEG$Name)

limma_GSE9750_DEG_rm=unique(limma_GSE9750_DEG$symbol)
SAM_GSE9750_DEG_rm=unique(SAM_GSE9750_DEG$symbol)
t_GSE9750_DEG_rm=unique(t_GSE9750_DEG$symbol)

GSE9750_DEG_shared=intersect(limma_GSE9750_DEG_rm,SAM_GSE9750_DEG_rm)
GSE9750_DEG_shared=intersect(GSE9750_DEG_shared,t_GSE9750_DEG_rm)

write.csv(limma_GSE9750_DEG,file = 'limma_GSE9750_DEG_symbol.csv')
write.csv(SAM_GSE9750_DEG,file = 'SAM_GSE9750_DEG_symbol.csv')
write.csv(t_GSE9750_DEG,file = 't_GSE9750_DEG_symbol.csv')
saveRDS(GSE9750_DEG_shared,file = 'GSE9750_DEG_shared.rds')
write.csv(GSE9750_DEG_shared, file = "GSE9750_DEG_shared.csv", row.names = FALSE)

# GSE168244 #################################

library(dplyr)
DESeq2_DEG = read.csv('DESeq2_GSE168244_DEG.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE168244_DEG.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE168244_DEG.csv', header = TRUE, row.names = 1)

DESeq2_DEG_rm=unique(rownames(DESeq2_DEG))
edgeR_DEG_rm=unique(rownames(edgeR_DEG))
NOISeq_DEG_rm=unique(rownames(NOISeq_DEG))

GSE168244_DEG_shared=intersect(DESeq2_DEG_rm,edgeR_DEG_rm)
GSE168244_DEG_shared=intersect(GSE168244_DEG_shared,NOISeq_DEG_rm)
saveRDS(GSE168244_DEG_shared,file = 'GSE168244_DEG_shared.rds')
write.csv(GSE168244_DEG_shared, file = "GSE168244_DEG_shared.csv", row.names = FALSE)


# GSE149450 ################################
#Bulk RNA:GSE149450 Gene ID conversion (https://www.biotools.fr/human/ensembl_symbol_converter)

library(dplyr)
DESeq2_DEG = read.csv('DESeq2_GSE149450_DEG_symbol.csv', header = TRUE,row.names=1)
edgeR_DEG = read.csv('edgeR_GSE149450_DEG_symbol.csv',  header = TRUE, row.names = 1)
NOISeq_DEG = read.csv('NOISeq_GSE149450_DEG_symbol.csv', header = TRUE, row.names = 1)

DESeq2_DEG=DESeq2_DEG %>% filter(symbol != '')
edgeR_DEG=edgeR_DEG %>% filter(symbol != '')
NOISeq_DEG=NOISeq_DEG %>% filter(symbol != '')

DESeq2_DEG_rm=unique(DESeq2_DEG$symbol)
edgeR_DEG_rm=unique(edgeR_DEG$symbol)
NOISeq_DEG_rm=unique(NOISeq_DEG$symbol)

GSE149450_DEG_shared=intersect(DESeq2_DEG_rm,edgeR_DEG_rm)
GSE149450_DEG_shared=intersect(GSE149450_DEG_shared,NOISeq_DEG_rm)
saveRDS(GSE149450_DEG_shared,file = 'GSE149450_DEG_shared.rds')
write.csv(GSE149450_DEG_shared, file = "GSE149450_DEG_shared.csv", row.names = FALSE)

# 4 datasets shared DEG ################################
DEG_shared1=intersect(GSE7803_DEG_shared,GSE9750_DEG_shared)
DEG_shared2=intersect(GSE168244_DEG_shared,GSE149450_DEG_shared)
DEG_shared=intersect(DEG_shared1,DEG_shared2)
write.csv(DEG_shared1,file = 'Microarray_DEG_shared.csv', row.names = FALSE)
write.csv(DEG_shared2,file = 'RNAseq_DEG_shared.csv', row.names = FALSE)
write.csv(DEG_shared,file = 'DEG_shared.csv', row.names = FALSE)



