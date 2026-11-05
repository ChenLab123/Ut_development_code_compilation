#figureS5a----------------------------------------------------------------------
library(pheatmap)
load('data/SMC_NMF50.Rdata')
programs.clusters <- data.frame(MP = geneNMF.metaprograms[["programs.clusters"]]) 
programsid <- row.names(programs.clusters)[programs.clusters$MP %in% c(1:18)]
jaccarddf <- geneNMF.metaprograms[["programs.jaccard"]][programsid,programsid]

library(RColorBrewer)
mp_colors <- c(
  '1'  = "#FDBF6F", 
  '2'  = "#E41A1C",  # 红
  '3'  = "#377EB8",  # 蓝
  '4'  = "#4DAF4A",  # 绿
  '5'  = "#984EA3",  # 紫
  '6'  = "#FF7F00",  # 橙
  '7'  = "#FFFF33",  # 黄
  '8'  = "#A65628",  # 棕
  '9'  = "#B3DE69", 
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

library(patchwork)
programs.anno <- data.frame(MP = as.character(programs.clusters[programs.clusters$MP %in% c(1:18),]))
rownames(programs.anno) <- rownames(programs.clusters)[programs.clusters$MP %in% c(1:18)]

p2 <- pheatmap(jaccarddf,cluster_rows = T,cluster_cols = T,annotation_legend = F,
               annotation_row = programs.anno ,annotation_colors = annotation_colors, 
               show_rownames = F,show_colnames = F,color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               cutree_rows = 10, cutree_cols = 10, treeheight_row = 0,breaks = seq(0, 0.6, length.out = 101),
               border_color = NA)

p2


pdf('figureS5/figS5a_MPcortotal_plot_label.pdf', width = 10, height = 6)
p2
dev.off()

#figureS5b----------------------------------------------------------------------
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


pdf('figureS5/figS5b_MPheatmap.pdf', width = 4, height = 6)
p1
dev.off()


#figureS5c----------------------------------------------------------------------
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


pdf('figureS5/figS5c_enrichment.pdf', width = 7, height = 5)
p1
dev.off()


#figureS5d----------------------------------------------------------------------
SMC_fa2 <- readRDS('data/SMC_fa2_pathway.RDS')
myeloid_fdg_coor <- SMC_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- SMC_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

pdf("figureS5/figS5d_SMC_fa2_pathway.pdf",width = 6, height = 6)
for (i in c("MAPKgs","PI3Kgs","Hippogs","Wntgs","ECMgs")) {
  p1<-ggplot(myeloid)+
    geom_point(aes_string(x = 'FDG1', y = 'FDG2',fill=i),shape = 21,stroke = 0.2,size = 3)+
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
    )+NoLegend()+
    labs(title = paste0('SMC_fa2_',i))
  print(p1)
}
dev.off()
#figureS5e----------------------------------------------------------------------
VSMC_fa2 <- readRDS('data/VSMC_fa2_pathway.RDS')
myeloid_fdg_coor <- VSMC_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- VSMC_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

pdf("figureS5/figS5e_VSMC_fa2_pathway.pdf",width = 6, height = 6)
for (i in c("MAPKgs","PI3Kgs","Hippogs","Wntgs","ECMgs")) {
  p1<-ggplot(myeloid)+
    geom_point(aes_string(x = 'FDG1', y = 'FDG2',fill=i),shape = 21,stroke = 0.2,size = 3)+
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
    )+NoLegend()+
    labs(title = paste0('SMC_fa2_',i))
  print(p1)
}
dev.off()

#figureS5f----------------------------------------------------------------------
library(reshape2)
library(ggplot2)

load('data/MPgenes.ucell.Rdata')
load('data/SMC_NMF50.Rdata')
table(colnames(combined.ucell) == colnames(seu))
seu$seurat_clusters <- combined.ucell$seurat_clusters

p1 <- DimPlot(seu,reduction = 'fa2',group.by = 'seurat_clusters',label = T)
ggsave(p1,filename = "figureS5/figS5f_fa2.pdf",width = 6, height = 6, units = "in", device='pdf')

#MP_Des_SMC
metadf <- seu@meta.data
metadf <- metadf[metadf$seurat_clusters %in% c('2','0','20','13','9','6'),]
metadf$seurat_clusters <- factor(metadf$seurat_clusters,levels = c('2','0','20','13','9','6'))
metadf$id <- row.names(metadf)

df <- melt(metadf,id.vars = c('id','seurat_clusters'),measure.vars = c("MP2","MP9","MP10","MP12"),
           value.name = 'ucell')

result <- aggregate(metadf, 
                    by=list(type=metadf$seurat_clusters), FUN=mean)

dfave <- melt(result,id.vars = c('type'),measure.vars = c("MP2","MP9","MP10","MP12"),
              value.name = 'ucell')

p1 <- ggplot(data = df, mapping = aes(x = seurat_clusters, y = ucell, colour = variable, fill = variable)) + 
  geom_violin(position="dodge", alpha=0.5,scale = "width")+theme_bw()+
  labs(x = 'DES+SMC lineage')+ 
  geom_line(data = dfave, mapping = aes(x = type, y = ucell, colour = variable, group = variable),size =1,position = position_dodge(width = 0.9))+ 
  geom_point(data = dfave, mapping = aes(x = type, y = ucell, colour = variable, group = variable),size =3,position = position_dodge(width = 0.9))

print(p1)
ggsave(p1,filename = 'figureS5/FigS5f_Des_SMC.pdf',width = 10, height = 6, units = "in", device='pdf')


#MP_ESC
metadf <- seu@meta.data
metadf <- metadf[metadf$seurat_clusters %in% c('2','0','20','13','16','15','4'),]
metadf$seurat_clusters <- factor(metadf$seurat_clusters,levels = c('2','0','20','13','16','15','4'))
metadf$id <- row.names(metadf)

df <- melt(metadf,id.vars = c('id','seurat_clusters'),measure.vars = c("MP2","MP9","MP10","MP12"),
           value.name = 'ucell')

result <- aggregate(metadf, 
                    by=list(type=metadf$seurat_clusters), FUN=mean)

dfave <- melt(result,id.vars = c('type'),measure.vars = c("MP2","MP9","MP10","MP12"),
              value.name = 'ucell')

p1 <- ggplot(data = df, mapping = aes(x = seurat_clusters, y = ucell, colour = variable, fill = variable)) + 
  geom_violin(position="dodge", alpha=0.5,scale = "width")+theme_bw()+
  labs(x = 'ESC lineage')+ 
  geom_line(data = dfave, mapping = aes(x = type, y = ucell, colour = variable, group = variable),size =1,position = position_dodge(width = 0.9))+ 
  geom_point(data = dfave, mapping = aes(x = type, y = ucell, colour = variable, group = variable),size =3,position = position_dodge(width = 0.9))

print(p1)
ggsave(p1,filename = 'figureS5/FigS5f_ESC.pdf',width = 10, height = 6, units = "in", device='pdf')


#MP_VSMC
metadf <- seu@meta.data
metadf <- metadf[metadf$seurat_clusters %in% c('1','8','5','12','14','3','11','7','10'),]
metadf$seurat_clusters <- factor(metadf$seurat_clusters,levels = c('1','8','5','12','14','3','11','7','10'))
metadf$id <- row.names(metadf)

df <- melt(metadf,id.vars = c('id','seurat_clusters'),measure.vars = c("MP2","MP9","MP10","MP12"),
           value.name = 'ucell')

result <- aggregate(metadf, 
                    by=list(type=metadf$seurat_clusters), FUN=mean)

dfave <- melt(result,id.vars = c('type'),measure.vars = c("MP2","MP9","MP10","MP12"),
              value.name = 'ucell')


p1 <- ggplot(data = df, mapping = aes(x = seurat_clusters, y = ucell, colour = variable, fill = variable)) + 
  geom_violin(position="dodge", alpha=0.5,scale = "width")+theme_bw()+
  labs(x = 'VSMC lineage')+ 
  geom_line(data = dfave, mapping = aes(x = type, y = ucell, colour = variable, group = variable),size =1,position = position_dodge(width = 0.9))+ 
  geom_point(data = dfave, mapping = aes(x = type, y = ucell, colour = variable, group = variable),size =3,position = position_dodge(width = 0.9))

print(p1)
ggsave(p1,filename = 'figureS5/FigS5f_VSMC.pdf',width = 12, height = 6, units = "in", device='pdf')

