####### GO- BP ###############
#DAVID
rm(list = ls())
gc()
library(tidyverse)
BP=read.delim('BP_DAVID.txt',header = T)
BP=BP %>% separate(Term,into = c('GO ID','Terms'),sep='~')
BP_DAVID=BP[,c('GO ID','Terms','PValue','Genes')]
write.csv(BP_DAVID,file = 'BP_DAVID_table.csv',row.names = F)

#enrichr
library(dplyr)
BP=read.delim('BP_enrichr.txt',header = T)
BP=BP %>% separate(Term,into = c('Terms','GO ID'),sep='\\(')
BP$`GO ID`= sub("\\)", "", BP$`GO ID`)
BP_enrichr=BP %>% filter(Adjusted.P.value<0.05)
BP_enrichr=BP_enrichr[,c('GO ID','Terms','Adjusted.P.value','Genes')]
colnames(BP_enrichr)=c('GO ID','Terms','PValue','Genes')
write.csv(BP_enrichr,file = 'BP_enrichr_table.csv',row.names = F)

#webgestalt
BP=read.delim('BP_WebGestalt.txt',header = T)
#BP=BP %>% filter(pValue<0.05)
BP_webgestalt=BP[,c('geneSet','description','pValue','userId')]#
colnames(BP_webgestalt)=c('GO ID','Terms','PValue','Genes')
write.csv(BP_webgestalt,file = 'BP_webgestalt_table.csv',row.names = F)


####### GO- CC ###############
#DAVID
rm(list = ls())
gc()
CC = read.delim('CC_DAVID.txt', header = T)
CC = CC %>% separate(Term, into = c('GO ID', 'Terms'), sep = '~')
CC = CC %>% filter(PValue<0.05)
CC_DAVID = CC[, c('GO ID', 'Terms', 'PValue', 'Genes')]
write.csv(CC_DAVID, file = 'CC_DAVID_table.csv', row.names = F)

#enrichr
library(dplyr)
CC = read.delim('CC_enrichr.txt', header = T)
CC = CC %>% separate(Term, into = c('Terms', 'GO ID'), sep = '\\(')
CC$`GO ID` = sub("\\)", "", CC$`GO ID`)
CC_enrichr = CC %>% filter(Adjusted.P.value < 0.05)
CC_enrichr = CC_enrichr[, c('GO ID', 'Terms', 'Adjusted.P.value', 'Genes')]
colnames(CC_enrichr) = c('GO ID', 'Terms', 'PValue', 'Genes')
write.csv(CC_enrichr, file = 'CC_enrichr_table.csv', row.names = F)

#webgestalt
CC = read.delim('CC_WebGestalt.txt', header = TRUE)
CC_webgestalt = CC[, c('geneSet', 'description', 'pValue', 'userId')]
colnames(CC_webgestalt) = c('GO ID', 'Terms', 'PValue', 'Genes')
write.csv(CC_webgestalt, file = 'CC_webgestalt_table.csv', row.names = FALSE)

####### GO- MF ###############
#DAVID
rm(list = ls())
gc()
MF = read.delim('MF_DAVID.txt', header = T)
MF = MF %>% separate(Term, into = c('GO ID', 'Terms'), sep = '~')
MF = MF %>% filter(PValue<0.05)
MF_DAVID = MF[, c('GO ID', 'Terms', 'PValue', 'Genes')]
write.csv(MF_DAVID, file = 'MF_DAVID_table.csv', row.names = F)

#enrichr
library(dplyr)
MF = read.delim('MF_enrichr.txt', header = T)
MF = MF %>% separate(Term, into = c('Terms', 'GO ID'), sep = '\\(')
MF$`GO ID` = sub("\\)", "", MF$`GO ID`)
MF_enrichr = MF %>% filter(Adjusted.P.value < 0.05)
MF_enrichr = MF_enrichr[, c('GO ID', 'Terms', 'Adjusted.P.value', 'Genes')]
colnames(MF_enrichr) = c('GO ID', 'Terms', 'PValue', 'Genes')
write.csv(MF_enrichr, file = 'MF_enrichr_table.csv', row.names = F)

# webgestalt
MF = read.delim('MF_WebGestalt.txt', header = TRUE)
#MF = MF %>% filter(pValue<0.05)
MF_webgestalt = MF[, c('geneSet', 'description', 'pValue', 'userId')]
colnames(MF_webgestalt) = c('GO ID', 'Terms', 'PValue', 'Genes')
write.csv(MF_webgestalt, file = 'MF_webgestalt_table.csv', row.names = FALSE)

####### KEGG ###############
#DAVID
rm(list = ls())
gc()
KEGG = read.delim('KEGG_DAVID.txt', header = T)
KEGG = KEGG %>% separate(Term, into = c('ID', 'Terms'), sep = ':')
KEGG = KEGG %>% filter(PValue < 0.05)
KEGG_DAVID = KEGG[, c('ID', 'Terms', 'PValue', 'Genes')]
write.csv(KEGG_DAVID, file = 'KEGG_DAVID_table.csv', row.names = F)

#enrichr
library(dplyr)
KEGG = read.delim('KEGG_enrichr.txt', header = T)
KEGG_enrichr = KEGG %>% filter(Adjusted.P.value < 0.05)
KEGG_enrichr = KEGG_enrichr[, c('Term', 'Adjusted.P.value', 'Genes')]
colnames(KEGG_enrichr) = c('Terms', 'PValue', 'Genes')
write.csv(KEGG_enrichr, file = 'KEGG_enrichr_table.csv', row.names = F)

# webgestalt
KEGG = read.delim('KEGG_WebGestalt.txt', header = TRUE)
KEGG_webgestalt = KEGG[, c('geneSet', 'description', 'pValue', 'userId')]
colnames(KEGG_webgestalt) = c('ID', 'Terms', 'PValue', 'Genes')
write.csv(KEGG_webgestalt, file = 'KEGG_webgestalt_table.csv', row.names = FALSE)

####### REACTOME ###############
#DAVID
rm(list = ls())
gc()
REACTOME = read.delim('REACTOME_DAVID.txt', header = T)
REACTOME = REACTOME %>% separate(Term, into = c('ID', 'Terms'), sep = '~')
REACTOME = REACTOME %>% filter(PValue < 0.05)
REACTOME_DAVID = REACTOME[, c('ID', 'Terms', 'PValue', 'Genes')]
write.csv(REACTOME_DAVID, file = 'REACTOME_DAVID_table.csv', row.names = F)

#enrichr
library(dplyr)
REACTOME = read.delim('REACTOME_enrichr.txt', header = T)
REACTOME = REACTOME %>% separate(Term, into = c('Terms', 'ID'), sep = 'R-')
REACTOME$ID = paste0('R-', REACTOME$ID)
REACTOME_enrichr = REACTOME %>% filter(Adjusted.P.value < 0.05)
REACTOME_enrichr = REACTOME_enrichr[, c('ID', 'Terms', 'Adjusted.P.value', 'Genes')]
colnames(REACTOME_enrichr) = c('ID', 'Terms', 'PValue', 'Genes')
write.csv(REACTOME_enrichr, file = 'REACTOME_enrichr_table.csv', row.names = F)

# webgestalt
REACTOME = read.delim('REACTOME_WebGestalt.txt', header = TRUE)
REACTOME_webgestalt = REACTOME[, c('geneSet', 'description', 'pValue', 'userId')]
colnames(REACTOME_webgestalt) = c('ID', 'Terms', 'PValue', 'Genes')
write.csv(REACTOME_webgestalt, file = 'REACTOME_webgestalt_table.csv', row.names = FALSE)


####### WikiPathways ###############
#DAVID
rm(list = ls())
gc()
WikiPathways = read.delim('WikiPathways_DAVID.txt', header = T)
WikiPathways = WikiPathways %>% separate(Term, into = c('ID', 'Terms'), sep = ':')
WikiPathways = WikiPathways %>% filter(PValue < 0.05)
WikiPathways_DAVID = WikiPathways[, c('ID', 'Terms', 'PValue', 'Genes')]
write.csv(WikiPathways_DAVID, file = 'WikiPathways_DAVID_table.csv', row.names = F)

#enrichr
library(dplyr)
WikiPathways = read.delim('WikiPathways_enrichr.txt', header = T)
WikiPathways = WikiPathways %>% separate(Term, into = c('Terms', 'ID'), sep = 'WP')
WikiPathways$ID = paste0('WP', WikiPathways$ID)
WikiPathways_enrichr = WikiPathways %>% filter(Adjusted.P.value < 0.05)
WikiPathways_enrichr = WikiPathways_enrichr[, c('ID', 'Terms', 'Adjusted.P.value', 'Genes')]
colnames(WikiPathways_enrichr) = c('ID', 'Terms', 'PValue', 'Genes')
write.csv(WikiPathways_enrichr, file = 'WikiPathways_enrichr_table.csv', row.names = F)

# webgestalt
WikiPathways = read.delim('WikiPathways_WebGestalt.txt', header = TRUE)
WikiPathways_webgestalt = WikiPathways[, c('geneSet', 'description', 'pValue', 'userId')]
colnames(WikiPathways_webgestalt) = c('ID', 'Terms', 'PValue', 'Genes')
write.csv(WikiPathways_webgestalt, file = 'WikiPathways_webgestalt_table.csv', row.names = FALSE)

