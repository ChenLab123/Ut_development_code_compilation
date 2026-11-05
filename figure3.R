library(Seurat)
library(ggplot2)

SMC_fa2 <-  readRDS('data/SMC_fa2.RDS')
#figure3a fa2---------------------------------------------------------------------

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

myeloid_fdg_coor <- SMC_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- SMC_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=subanno),shape = 21,stroke = 0.2,size = 3)+
  theme_bw()+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )+NoLegend()
p1
ggsave(p1,filename = "figure3/fig3a_SMC_fa2.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure3b cytotrace---------------------------------------------------------------------
load('data/SMC_CytoTRACE.Rdata')
cytodf <- as.data.frame(results[["CytoTRACE"]])
SMC_cyto <- AddMetaData(SMC_fa2,cytodf,col.name = "CytoTrace")

myeloid_fdg_coor <- SMC_cyto@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- SMC_cyto@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=CytoTrace),shape = 21,stroke = 0.2,size = 3)+
  theme_bw()+
  scale_fill_continuous(palette = c("#F2F9FD","#E5F0F9", "#C4E2F3", "#90CAE9", "#55AEDB"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )+NoLegend()
p1
ggsave(p1,filename = "figure3/fig3b_cyto.pdf",width = 6, height = 6, units = "in", device='pdf')

p1<-ggplot(myeloid)+
  geom_boxplot(aes(x = reorder(subanno, -CytoTrace), y = CytoTrace,fill=subanno))+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme_classic()+
  NoLegend()+
  theme(axis.title.x=element_blank())

p1
ggsave(p1,filename = "figure3/fig3b_cyto_box.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure3b scVelocity---------------------------------------------------------------------
library(anndata)
library(reticulate)
ad <- import("anndata")
final_ad <- ad$read_h5ad("data/ut1-27_velocity.h5ad")
vedf <- final_ad$obs['velocity_pseudotime']

myeloid_fdg_coor <- read.csv('data/cell_embeddings.csv',row.names = 1)
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- read.csv('data/cell_clusters.csv',row.names = 1)
myeloid <- merge(myeloid_fdg_coor,vedf,by = 0)
myeloid <- merge(myeloid,metadf,by.x = 'Row.names',by.y = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=velocity_pseudotime),shape = 21,stroke = 0.2,size = 3)+
  theme_bw()+
  scale_fill_continuous(palette = c("#F2F9FD","#E5F0F9", "#C4E2F3", "#90CAE9", "#55AEDB"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )+NoLegend()
p1
ggsave(p1,filename = "figure3/fig3b_scve.pdf",width = 6, height = 6, units = "in", device='pdf')

p1<-ggplot(myeloid)+
  geom_boxplot(aes(x = reorder(subanno, velocity_pseudotime), y = velocity_pseudotime,fill=subanno))+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme_classic()+
  NoLegend()+
  theme(axis.title.x=element_blank())

p1

ggsave(p1,filename = "figure3/fig3b_scve_box.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure3b monocle3---------------------------------------------------------------------
load('data/SMC_monocle3.Rdata')
SMC_mono <- AddMetaData(
  object = SMC_fa2,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

myeloid_fdg_coor <- SMC_mono@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- SMC_mono@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=pseudotime),shape = 21,stroke = 0.2,size = 3)+
  theme_bw()+
  scale_fill_continuous(palette = c("#F2F9FD","#E5F0F9", "#C4E2F3", "#90CAE9", "#55AEDB"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )+NoLegend()
p1
ggsave(p1,filename = "figure3/fig3b_mono3.pdf",width = 6, height = 6, units = "in", device='pdf')

p1<-ggplot(myeloid)+
  geom_boxplot(aes(x = reorder(subanno, pseudotime), y = pseudotime,fill=subanno))+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme_classic()+
  NoLegend()+
  theme(axis.title.x=element_blank())

p1

ggsave(p1,filename = "figure3/fig3b_mono3_box.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure3b slingshot---------------------------------------------------------------------
load('data/SMC_slingshot.Rdata')
slingpe <- sling_obj@assays@data@listData[["pseudotime"]]
SMC_sling <- AddMetaData(
  object = SMC_fa2,
  metadata = slingpe,
  col.name = c("slingshot1","slingshot2")
)

myeloid_fdg_coor <- SMC_sling@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- SMC_sling@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

myeloid$slingshot <- ifelse(is.na(myeloid$slingshot2),myeloid$slingshot1,myeloid$slingshot2)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=slingshot),shape = 21,stroke = 0.2,size = 3)+
  theme_bw()+
  scale_fill_continuous(palette = c("#F2F9FD","#E5F0F9", "#C4E2F3", "#90CAE9", "#55AEDB"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )+NoLegend()
p1

ggsave(p1,filename = "figure3/fig3b_slingshot.pdf",width = 6, height = 6, units = "in", device='pdf')

p1<-ggplot(myeloid)+
  geom_boxplot(aes(x = reorder(subanno, slingshot), y = slingshot,fill=subanno))+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme_classic()+
  NoLegend()+
  theme(axis.title.x=element_blank())

p1

ggsave(p1,filename = "figure3/fig3b_slingshot_box.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure3c cor---------------------------------------------------------------------
library(corrplot)
scRNA <- readRDS('data/SMC_fa2.RDS')
scRNAid <- scRNA@meta.data
scRNAid$orig.ident <-  gsub('Ut','ut_',scRNAid$orig.ident )
scRNAid$cellid <- gsub('Ut26_|Ut13_|Ut10.1_|Ut10.2_|Ut21_|Ut23_|Ut27_','',row.names(scRNAid))
scRNAid$cellid <- gsub('-.*','x',scRNAid$cellid)
scRNAid$cellid <- paste0(scRNAid$orig.ident,':',scRNAid$cellid)
scRNAid$cellid <- gsub('ut_10.1','ut_10_1',scRNAid$cellid)
scRNAid$cellid <- gsub('ut_10.2','ut_10_2',scRNAid$cellid)
scRNAid$cellid <- gsub('ut_27','Ut_27-3RNA',scRNAid$cellid)
scRNAid <- scRNAid[,c("subanno","cellid" )]

library(anndata)
library(reticulate)
ad <- import("anndata")
final_ad <- ad$read_h5ad("data/ut1-27_velocity.h5ad")
vedf <- final_ad$obs['velocity_pseudotime']

load('data/SMC_CytoTRACE.Rdata')
cytodf <- as.data.frame(results[["CytoTRACE"]])

load('data/SMC_monocle3.Rdata')
mono3df <- as.data.frame(cds@principal_graph_aux@listData$UMAP$pseudotime)

load('data/SMC_slingshot.Rdata')
slingpe <- sling_obj@assays@data@listData[["pseudotime"]]

scRNAid <- merge(scRNAid,cytodf,by = 0)
scRNAid <- merge(scRNAid,mono3df,by.x = 'Row.names',by.y =0)
scRNAid <- merge(scRNAid,slingpe,by.x = 'Row.names',by.y =0)
scRNAid <- merge(scRNAid,vedf,by.x ='cellid',by.y = 0,all.x = T)


colnames(scRNAid)[4:8] <- c('CytoTRACE','monocle3','slingshot1','slingshot2','scVelocity')
write.csv(scRNAid,'data/SMC_pe.csv',row.names = F)

scRNAid <- read.csv('data/SMC_pe.csv')
scRNAid <- scRNAid[!is.na(scRNAid$slingshot1),]
cordf <- scRNAid[,c(4,5,6,8)]

cor_matrix <- cor(cordf,use = "pairwise.complete.obs")

pdf('figure3/fig3c_cor_slingshot1.pdf', width = 6, height = 6)
corrplot.mixed(cor_matrix, 
               upper = "circle", 
               lower = "number",
               tl.col = "black")
dev.off()

#figure3d----------------------------------------------------------------------
library(splines)
library(factoextra)
library(FactoMineR)
library(scales)
library(dplyr)

#slingshot1 DES SMC--------------------------------------------------------------------------
SMCpe <- read.csv('data/SMC_pe.csv')
SMCpe <- SMCpe[!is.na(SMCpe$slingshot1),]
SMCpe <- arrange(SMCpe, slingshot1)

load(file = 'data/exp500_slingshot1.Rdata')

order1 = arrange(exp500, atan2.1)

library(ComplexHeatmap)
library(RColorBrewer)
library(grid)

SMCpedf <- SMCpe
SMCpedf$CytoTRACE <- rescale(SMCpedf$CytoTRACE, to = c(0, 10))
SMCpedf$monocle3 <- rescale(SMCpedf$monocle3, to = c(0, 10))
SMCpedf$scVelocity <- rescale(SMCpedf$scVelocity, to = c(0, 10))
SMCpedf$slingshot1 <- rescale(SMCpedf$slingshot1, to = c(0, 10))
SMCpedf <- SMCpedf[,c('CytoTRACE','monocle3','scVelocity','slingshot1')] 
SMCpedf <- t(SMCpedf)

p1 = pheatmap(order1[,1:500], cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = FALSE, show_colnames = FALSE,
              col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
              breaks = seq(-1.5, 1.5, length.out = 101))
p1

p2 = pheatmap(SMCpedf,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,
              color = colorRampPalette(colors = c("white", "#6f8fc6"))(100), gaps_row = c(1,2,3))

p2


pdf('figure3/fig3d_MPgenes_DES.SMC.pdf', width = 8, height = 6)
p1
dev.off()


pdf('figure3/fig3d_petime_DES.SMC.pdf', width = 8, height = 1)
p2
dev.off()

write.csv(order1,file = 'data/MPgenes_DES.SMC.csv')


#slingshot2 ESC----------------------------------------------------------------------
SMCpe <- read.csv('data/SMC_pe.csv')
SMCpe <- SMCpe[!is.na(SMCpe$slingshot2),]
SMCpe <- arrange(SMCpe, slingshot2)
load(file = 'data/exp500_slingshot2.Rdata')

# 按角度值对基因进行排序，获得4种顺序
order1 = arrange(exp500, atan2.1)

library(ComplexHeatmap)
library(RColorBrewer)
library(grid)

SMCpedf <- SMCpe
SMCpedf$CytoTRACE <- rescale(SMCpedf$CytoTRACE, to = c(0, 10))
SMCpedf$monocle3 <- rescale(SMCpedf$monocle3, to = c(0, 10))
SMCpedf$scVelocity <- rescale(SMCpedf$scVelocity, to = c(0, 10))
SMCpedf$slingshot2 <- rescale(SMCpedf$slingshot2, to = c(0, 10))
SMCpedf <- SMCpedf[,c('CytoTRACE','monocle3','scVelocity','slingshot2')] 
SMCpedf <- t(SMCpedf)

p1 = pheatmap(order1[,1:500], cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = FALSE, show_colnames = FALSE,
              col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
              breaks = seq(-1.5, 1.5, length.out = 101))
p1

p2 = pheatmap(SMCpedf,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,
              color = colorRampPalette(colors = c("white", "#FFAA00"))(100), gaps_row = c(1,2,3))

p2

pdf('figure3/fig3d_MPgenes_ESC.pdf', width = 8, height = 6)
p1
dev.off()

pdf('figure3/fig3d_petime_ESC.pdf', width = 8, height = 1)
p2
dev.off()

write.csv(order1,file = 'data/MPgenes_ESC.csv')

#MP plot--------------------------------------------------------------------------
library(UCell)
library(Seurat)
library(pheatmap)
library(scales)
library(dplyr)

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

load('data/SMC_NMF50.Rdata')

load('data/SMC_slingshot.Rdata')
slingpe <- sling_obj@assays@data@listData[["pseudotime"]]
SMC_sling <- AddMetaData(
  object = seu,
  metadata = slingpe,
  col.name = c("slingshot1","slingshot2")
)


metadf <- SMC_sling@meta.data

metadf <- metadf[!is.na(metadf$slingshot1),]
metadf <- arrange(metadf, slingshot1)

metadf$MP10 <- rescale(metadf$MP10, to = c(0, 10))
metadf$MP12 <- rescale(metadf$MP12, to = c(0, 10))
metadf$slingshot1 <- rescale(metadf$slingshot1, to = c(0, 10))
metadf1 <- metadf[,c('slingshot1','subanno')] 
metadf <- metadf[,c('MP10','MP12','slingshot1')] 
metadf <- t(metadf)

p3 = pheatmap(metadf,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,
              color = colorRampPalette(colors = c("white", "#6f8fc6"))(100), gaps_row = c(1,2))

p3


pdf('figure3/fig3d_MPbar_DES.SMC.pdf', width = 8, height = 1)
p3
dev.off()


p4 = ggplot(metadf1, aes(x = slingshot1, y = 1, fill = subanno)) +
  geom_tile(aes(width = 0.9, height = 0.6), color = NA) +
  scale_fill_manual(values = SEMPalette) +
  theme_void() +
  theme(legend.position = "bottom")


pdf('figure3/fig3d_MPbar_DES.SMC_1.pdf', width = 8, height = 0.5)
print(p4) 
dev.off()


metadf <- SMC_sling@meta.data
metadf <- metadf[!is.na(metadf$slingshot2),]
metadf <- arrange(metadf, slingshot2)

metadf$MP10 <- rescale(metadf$MP10, to = c(0, 10))
metadf$MP9 <- rescale(metadf$MP9, to = c(0, 10))
metadf$slingshot2 <- rescale(metadf$slingshot2, to = c(0, 10))
metadf1 <- metadf[,c('slingshot2','subanno')]
metadf <- metadf[,c('MP10','MP9','slingshot2')] 
metadf <- t(metadf)

p3 = pheatmap(metadf,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,
              color = colorRampPalette(colors = c("white", "#FFAA00"))(100), gaps_row = c(1,2))

p3


pdf('figure3/fig3d_MPbar_ESC.pdf', width = 8, height = 1.5)
p3
dev.off()

p4 =  ggplot(metadf1, aes(x = slingshot2, y = 1, fill = subanno)) +
  geom_tile(aes(width = 0.9, height = 0.6), color = NA) +
  scale_fill_manual(values = SEMPalette) +
  theme_void() +
  theme(legend.position = "bottom")

pdf('figure3/fig3d_MPbar_ESC_1.pdf', width = 8, height = 0.5)
print(p4) 
dev.off()


#figure3g heatmap ----------------------------------------------------------------------
library(pheatmap)
load('data/SMC_NMF50.Rdata')
programs.clusters <- data.frame(MP = geneNMF.metaprograms[["programs.clusters"]]) 
programsid <- row.names(programs.clusters)[programs.clusters$MP %in% c(4,10,12,13,18)]
jaccarddf <- geneNMF.metaprograms[["programs.jaccard"]][programsid,programsid]

library(RColorBrewer)
mp_colors <- c(
  '2'  = "#E41A1C",  # 红
  '3'  = "#377EB8",  # 蓝
  '4'  = "#4DAF4A",  # 绿
  '5'  = "#984EA3",  # 紫
  '6'  = "#FF7F00",  # 橙
  '7'  = "#FFFF33",  # 黄
  '8'  = "#A65628",  # 棕
  '10' = "#F781BF",  # 粉
  '11' = "#999999",  # 灰
  '12' = "#66C2A5",  # 青绿
  '13' = "#FC8D62",  # 橙红
  '14' = "#8DA0CB",  # 浅蓝
  '15' = "#E78AC3",  # 浅粉
  '16' = "#A6D854",  # 草绿
  '17' = "#FFD92F",  # 金黄
  '18' = "#E5C494"   # 米色
)
annotation_colors <- list(MP = mp_colors)

p1 <- pheatmap(jaccarddf,cluster_rows = T,cluster_cols = T,
               show_rownames = F,show_colnames = F,color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               cutree_rows = 5, cutree_cols = 5, treeheight_row = 0,breaks = seq(0, 0.6, length.out = 101),
               border_color = NA)
p1

library(patchwork)
programs.anno <- data.frame(MP = as.character(programs.clusters[programs.clusters$MP %in% c(4,10,12,13,18),]))
rownames(programs.anno) <- rownames(programs.clusters)[programs.clusters$MP %in% c(4,10,12,13,18)]


p2 <- pheatmap(jaccarddf,cluster_rows = T,cluster_cols = T,annotation_legend = F,
               annotation_row = programs.anno ,annotation_colors = annotation_colors, 
               show_rownames = F,show_colnames = F,color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               cutree_rows = 5, cutree_cols = 5, treeheight_row = 0,breaks = seq(0, 0.6, length.out = 101),
               border_color = NA)

p2


pdf('figure3/fig3g_MPcorsub_plot_label.pdf', width = 10, height = 6)
p2
dev.off()



#figure3h----------------------------------------------------------------------
library(pheatmap)
library(tidyr)
library(Seurat)
library(UCell)
library(dplyr)
library(reshape2)

SEMPalette =  c("#ffadad",'#2c6e49', '#fefee3', '#ffc9b9', '#4c956c', '#d68c45') 
names(SEMPalette) <- c("ESC", "DES.SMC",'Mesen0','Mesen6',"VSMC",'Mesen2')

load('data/SMC_NMF50.Rdata')
mp.genes <- geneNMF.metaprograms$metaprograms.genes

metadf <- seu@meta.data
metadf$id <- row.names(metadf)
metadf <- metadf[,c('subanno',names(mp.genes))]

result <- aggregate(metadf, 
                    by=list(metadf$subanno), FUN=mean)
row.names(result) <- result$Group.1
result <- result[,-1]
result <- result[,-1]
result <- scale(result)
result <- t(result)
result <- result[,c('Mesen0','Mesen6', "DES.SMC","ESC",'Mesen2',"VSMC")]

annotation_df <- data.frame(
  Group = factor(c('Mesen0','Mesen6', "DES.SMC","ESC",'Mesen2',"VSMC"),
                 levels = c('Mesen0','Mesen6', "DES.SMC","ESC",'Mesen2',"VSMC")# 显式顺序
  )
)
rownames(annotation_df) <- c('Mesen0','Mesen6', "DES.SMC","ESC",'Mesen2',"VSMC")

# 自定义颜色映射
ann_colors <- list(
  Group = SEMPalette
)


p1 <- pheatmap(result,cluster_rows = F,cluster_cols = F,
               show_rownames = T,show_colnames = F,
               color = colorRampPalette(colors = c("#6f8fc6", "white", "#ff6361"))(100),
               treeheight_row = 0,treeheight_col = 0,
               annotation_col = annotation_df,annotation_colors = ann_colors,
               annotation_names_col = T,
               breaks = seq(-1.5, 1.5, length.out = 101))
p1


pdf('figure3/fig3h_MPheatmap.pdf', width = 4, height = 6)
p1
dev.off()


#figure3i----------------------------------------------------------------------
library(pheatmap)
library(tidyr)
kegglist <- read.csv('data/MP_KEGG1_NMF50.csv')
kegglist$MP <- paste0('MP',kegglist$MP)
kegg_width <- spread(kegglist, key = "MP",
                     value = "pvalue")

row.names(kegg_width) <- kegg_width$Description
kegg_width <- kegg_width[,-1]
kegg_width <- kegg_width[unique(kegglist$Description),unique(kegglist$MP)]

kegg_width[is.na(kegg_width)] <- 1
kegg_width <- -log10(kegg_width)


p1 <- pheatmap(kegg_width,cluster_rows = F,cluster_cols = F,
               show_rownames = T,show_colnames = T,
               color = colorRampPalette(colors = c( "white", "red"))(100),
               treeheight_row = 0,treeheight_col = 0,
               breaks = seq(0, 1.5, length.out = 101))
p1


pdf('figure3/fig3i_enrichment.pdf', width = 7, height = 5)
p1
dev.off()

