#### top 10 of shared 22 ####
rm(list = ls())
gc()
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE149450")
Data1=readRDS('GSE149450_enrich.rds')
top10=Data1[order(abs(Data1$mean_nes),decreasing = T),]
top10=top10[1:10,]
GSE149450=readRDS('GSE149450_NES_2.rds')
ov=intersect(GSE149450$pathway,top10$pathway)

#### top 10 of shared 203  ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE168244")
Data1=readRDS('GSE168244_enrich.rds')
top10=Data1[order(abs(Data1$mean_nes),decreasing = T),]
top10=top10[1:10,]
GSE149450=readRDS('GSE168244_NES_2.rds')
ove=intersect(GSE149450$pathway,top10$pathway)


#### top 10 of shared 735 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE9750")
Data1=readRDS('GSE9750_enrich.rds')
top10=Data1[order(abs(Data1$mean_nes),decreasing = T),]
top10=top10[1:10,]
GSE149450=readRDS('GSE9750_NES_2.rds')
over=intersect(GSE149450$pathway,top10$pathway)


#### top 10 of shared 596 ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/GSE7803")
Data1=readRDS('GSE7803_enrich.rds')
top10=Data1[order(abs(Data1$mean_nes),decreasing = T),]
top10=top10[1:10,]
GSE149450=readRDS('GSE7803_NES_2.rds')
overl=intersect(GSE149450$pathway,top10$pathway)
overla=c()
overla=c(ov,ove,over,overl)
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
saveRDS(overla,file = 'Top10_38_pathway_name.rds')

#### Sankey data ####
setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
overla=readRDS('Top10_38_pathway_name.rds')

result_NES_2 <- readRDS("~/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/result_NES_2.rds")
result_NES_2=result_NES_2[result_NES_2$Pathway%in% overla,]
Sankey=result_NES_2[,c(1,9,10)]
KG=c('CDK1','TOP2A','BRCA1','CDC20','FEN1','BUB1B','CDC45','KIF23','KIF11','MCM7')

filterpathway=c()
for (i in 1:nrow(Sankey)){
  show=intersect(KG,unlist(result_NES_2[i, 8]))
  if(length(show)>0){
    add=i
    if(i==1){
      filterpathway=add
    }else{
      filterpathway=c(filterpathway,add)
    }}
}

Sankey=result_NES_2[filterpathway,]
Sankey=Sankey[which(Sankey$NES>2.22),]
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
#saveRDS(Sankey_data,file = 'Sankey_data_NES_2_KG.rds')

Tb=as.data.frame(table(Sankey_data$Sample))
sankey.ed=Sankey_data[,c('Sample','Pathway','KeyGene')]
#write.csv(Sankey_data[,c('Sample','Pathway','KeyGene')],file = 'Sankey_data_NES_2_KG.csv',row.names = F)
Tb2=as.data.frame(table(Sankey_data$Pathway))
sankey.ed=sankey.ed[!(sankey.ed$Pathway %in% c("GABRIELY_MIR21_TARGETS", 
                                               "SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_UP",
                                               'DAZARD_RESPONSE_TO_UV_NHEK_DN',
                                               'REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS',
                                               'REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX')), ]




################ sankey plot #####
library(ggalluvial)
library(ggplot2)

workdir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea'
savedir='/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/Sankey plot'
setwd(workdir)

#读取数据
#sankey.ed <- read.csv('Sankey_data_NES_3_09_KG.csv', header = T)
# 绘图
sankey_colors<-c("#0072B2", "#999999","#D55E00", "#E69F00",'#cfe4ff',"#009E73", '#6f6f6f', "#56B4E9", '#b395bd', "#e9c46a")

p.sankey.ED <- ggplot(data = sankey.ed,
                      aes(axis1 = sankey.ed$Sample,
                          axis2 = sankey.ed$Pathway,
                          axis3 = sankey.ed$KeyGene)) +
  scale_x_discrete(limits = c("Sample", 'Pathway', "KeyGene")) +
  geom_alluvium(aes(fill = sankey.ed$KeyGene, alpha =1),
                curve_type = "arctangent") +
  scale_color_manual(values = sankey_colors)+
  scale_fill_manual(values = sankey_colors)+
  geom_stratum(alpha = 0,color = adjustcolor( "white", alpha.f = 1),size=1.2, fill = 'white') + 
  geom_text(stat = "stratum",cex=12, aes(label = after_stat(stratum))) +
  theme_void()+
  theme(legend.position="none",
        axis.text = element_text(size = 20),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid=element_blank())
p.sankey.ED


setwd(savedir)
ggsave(p.sankey.ED, filename = "Sankey.top10inter22-m-h34-w27.pdf", width = 27, height = 34)














##### NES>3 sankey data ##################################
# setwd("/Users/wuhongbo/Desktop/数据分析/Cervical Cancer/middle_file/fgsea")
# result_NES_30 <- readRDS("~/Desktop/数据分析/Cervical Cancer/middle_file/fgsea/result_NES_3.rds")
# KG=c('CDK1','TOP2A','BRCA1','CDC20','FEN1','BUB1B','CDC45','KIF23','KIF11','MCM7')
# Sankey=result_NES_30[,c(1,9,10)]
# 
# filterpathway <- c()
# for (i in 1:nrow(Sankey)){
#   show=intersect(KG,unlist(result_NES_30[i, 8]))
#   if(length(show)>0){
#     add=i
#     if(i==1){
#       filterpathway=add
#     }else{
#       filterpathway=c(filterpathway,add)
#     }}
# }
# 
# Sankey=result_NES_30[filterpathway,]
# Sankey_data=data.frame()
# for (i in 1:nrow(Sankey)){
#   show=intersect(KG,unlist(Sankey[i, 8]))
#   new_rows=do.call(rbind,replicate(length(show),Sankey[i,],simplify = F))
#   new_rows$KeyGene=show
#   if(i==1){
#     Sankey_data=new_rows
#   }else{
#     Sankey_data=rbind(Sankey_data,new_rows)
#   }
# }
# colnames(Sankey_data)[1]='Pathway'
# #saveRDS(Sankey_data,file = 'Sankey_data_NES_3_0_KG.rds')
# #write.csv(Sankey_data[,c('Sample','Method','Pathway','KeyGene')],file = 'Sankey_data_NES_3_09_KG.csv',row.names = F)
# 
# 
# 
# 
# 
# 
# 









