#figureS2a SEM-----------------------------------------------------------------
SEM <- readRDS('data/SEM_subanno.RDS')
SEM$subanno <- ifelse(SEM$subanno %in% c('DES.SMC','VSMC'),SEM$subanno,as.character(SEM$harmony_anno27))
SEM <- subset(SEM ,subset = subanno == 'SMC',invert = T)

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

myeloid_fdg_coor <- SEM@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- SEM@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=subanno),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme_bw()+
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
ggsave(p1,filename = "figureS2/figS2a_SEM_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS2b Immu-----------------------------------------------------------------
Immu <- readRDS('data/Immu.RDS')
ImmuPalette =  c(
  "#b83c4f",  "#d96a78", "#f2a0ad", "#f5d0d8","#b6b1cf", "#8f8bc2",  "#6f8fc6", 
  "#4f77b8",  "#3e5aa8",  "#2f3e8f" ) 
names(ImmuPalette) <- c('T','B','NK','Fibrocyte','Progenitor','Macro.Eryt',
                        'Macro.monocyte-derived','Macro.TRM','Macro.prolif','Neutrophile')

myeloid_fdg_coor <- Immu@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- Immu@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=Immu_anno),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(ImmuPalette),breaks = names(ImmuPalette))+
  theme_bw()+
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
ggsave(p1,filename = "figureS2/figS2b_Immu_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS2c endo-----------------------------------------------------------------
endo <- readRDS('data/endo.RDS')
endoPalette =  c('#ef7c8e', '#fae8e0', '#b6e2d3', '#d8a7b1' ) 
names(endoPalette) <- c('AE','CE','Mural.Cell','VE')

myeloid_fdg_coor <- endo@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- endo@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=endo_anno),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(endoPalette),breaks = names(endoPalette))+
  theme_bw()+
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
ggsave(p1,filename = "figureS2/figS2c_endo_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS2d SEM-----------------------------------------------------------------
library(pheatmap)
SEM <- readRDS('data/SEM_subanno.RDS')
SEM$subanno <- ifelse(SEM$subanno %in% c('DES.SMC','VSMC'),SEM$subanno,as.character(SEM$harmony_anno27))
SEM <- subset(SEM ,subset = subanno == 'SMC',invert = T)

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

Idents(SEM) <- 'subanno'
use.gene <- c('Dpt', 'Igfbp6', 'Aebp1', 'Pdgfra','Gas5', 'Dlk1', 'Mfap4', 'Tgfbi',
              "Myh11", "Cnn1","Pcp4","Acta2","Gm13889","Mustn1","Rgs5" ,"Notch3",
              'ESR1', 'THY1', 'DES', 'COL1A1', 'COL1A2', 'LUM', 'DCN','FN1','DPP6')

plot.gene <- toupper(use.gene)

gene_exp <- AverageExpression(SEM,features = plot.gene,layer = 'data' )
gene_exp <- gene_exp[["RNA"]]
annotation_col <- data.frame(SEM_anno = gene_exp@Dimnames[[2]])
row.names(annotation_col) <- gene_exp@Dimnames[[2]]
ann_colors = list(SEM_anno = SEMPalette )


p1 <- pheatmap(gene_exp,  #要绘制热图的矩阵
               color = colorRampPalette(c('blue','white','red'))(100),
               border_color = "white", 
               scale = "row",
               cluster_rows = T,
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               fontsize = 8 ,
               annotation_legend = FALSE,
               annotation_col = annotation_col,
               annotation_colors = ann_colors
)

p1
ggsave(p1,filename = "figureS2/figS2d_SEM_heatmap.pdf",width = 4, height = 6, units = "in", device='pdf')


#figureS2e Immu-----------------------------------------------------------------
library(pheatmap)
Immu <- readRDS('data/Immu.RDS')

ImmuPalette =  c(
  "#b83c4f",  "#d96a78", "#f2a0ad", "#f5d0d8","#b6b1cf", "#8f8bc2",  "#6f8fc6", 
  "#4f77b8",  "#3e5aa8",  "#2f3e8f") 
names(ImmuPalette) <- c('T','B','NK','Fibrocyte','Progenitor','Macro.Eryt',
                        'Macro.monocyte-derived','Macro.TRM','Macro.prolif','Neutrophile')

Idents(Immu) <- 'Immu_anno'
Immumarker <- FindAllMarkers(Immu,min.pct = 0.3,logfc.threshold = 1,only.pos = T)
Immumarker <- Immumarker[Immumarker$p_val_adj < 0.05,]
use.gene <- Immumarker %>% group_by(cluster) %>% do(head(., n = 3))

plot.gene <- use.gene$gene

gene_exp <- AverageExpression(Immu,features = plot.gene,layer = 'data' )
gene_exp <- gene_exp[["RNA"]]
annotation_col <- data.frame(Immu_anno = gene_exp@Dimnames[[2]])
row.names(annotation_col) <- gene_exp@Dimnames[[2]]
ann_colors = list(Immu_anno = ImmuPalette )


p1 <- pheatmap(gene_exp,  #要绘制热图的矩阵
               color = colorRampPalette(c('blue','white','red'))(100),
               border_color = "white", 
               scale = "row",
               cluster_rows = T,
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               fontsize = 8 ,
               annotation_legend = FALSE,
               annotation_col = annotation_col,
               annotation_colors = ann_colors
)

p1
ggsave(p1,filename = "figureS2/figS2e_Immu_heatmap.pdf",width = 6, height = 10, units = "in", device='pdf')

#figureS2f endo-----------------------------------------------------------------
library(pheatmap)
endo <- readRDS('data/endo.RDS')

endoPalette =  c('#ef7c8e', '#fae8e0', '#b6e2d3', '#d8a7b1' ) 
names(endoPalette) <- c('AE','CE','Mural.Cell','VE')

Idents(endo) <- 'endo_anno'
endomarker <- FindAllMarkers(endo,min.pct = 0.3,logfc.threshold = 1,only.pos = T)
endomarker <- endomarker[endomarker$p_val_adj < 0.05,]
use.gene <- endomarker %>% group_by(cluster) %>% do(head(., n = 3))

plot.gene <- use.gene$gene

gene_exp <- AverageExpression(endo,features = plot.gene,layer = 'data' )
gene_exp <- gene_exp[["RNA"]]
annotation_col <- data.frame(endo_anno = gene_exp@Dimnames[[2]])
row.names(annotation_col) <- gene_exp@Dimnames[[2]]
ann_colors = list(endo_anno = endoPalette )


p1 <- pheatmap(gene_exp,  #要绘制热图的矩阵
               color = colorRampPalette(c('blue','white','red'))(100),
               border_color = "white", 
               scale = "row",
               cluster_rows = T,
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               fontsize = 8 ,
               annotation_legend = FALSE,
               annotation_col = annotation_col,
               annotation_colors = ann_colors
)

p1
ggsave(p1,filename = "figureS2/figS2f_endo_heatmap.pdf",width = 3, height = 5, units = "in", device='pdf')


#figureS2g SEM-----------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

SEM <- readRDS('data/SEM_subanno.RDS')
SEM$subanno <- ifelse(SEM$subanno %in% c('DES.SMC','VSMC'),SEM$subanno,as.character(SEM$harmony_anno27))
SEM <- subset(SEM ,subset = subanno == 'SMC',invert = T)

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')
tmp.meta = SEM@meta.data

tmp = tmp.meta[,c("time","subanno")]
tmp$new = paste(tmp$time,tmp$subanno,sep = "_")
rownames(tmp) = NULL
tmp2 = tmp[!duplicated(tmp),]
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(tmp[tmp$time == tmp2[t,"time"] & tmp$subanno == tmp2[t,"subanno"],])/nrow(tmp[tmp$time == tmp2[t,"time"],])
}

position <- read.csv('data/fig1d_alluvial_plot.csv')

tmp2 <- merge(tmp2,position,by = 'time', all =T)


p1 <- ggplot(tmp2,aes(x = position,stratum = subanno,alluvium = subanno,
                      y = rate,fill = subanno,label = subanno)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette)) +
  geom_flow(alpha = 0.3) +
  geom_stratum(width = 0.4,colour = "white", alpha = 0.7) +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        axis.text = element_text(size = 12),
        legend.position = "none")+ 
  scale_x_continuous(breaks = position$position,labels = position$time)
p1


pdf('figureS2/figS2g_alluvial_SEM.pdf', width = 10, height = 6)
p1
dev.off()

#figureS2h Immu-----------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

Immu <- readRDS('data/Immu.RDS')

ImmuPalette =  c(
  "#b83c4f",  "#d96a78", "#f2a0ad", "#f5d0d8","#b6b1cf", "#8f8bc2",  "#6f8fc6", 
  "#4f77b8",  "#3e5aa8",  "#2f3e8f") 
names(ImmuPalette) <- c('T','B','NK','Fibrocyte','Progenitor','Macro.Eryt',
                        'Macro.monocyte-derived','Macro.TRM','Macro.prolif','Neutrophile')

tmp.meta = Immu@meta.data

tmp = tmp.meta[,c("time","Immu_anno")]
tmp$new = paste(tmp$time,tmp$Immu_anno,sep = "_")
rownames(tmp) = NULL
tmp2 = tmp[!duplicated(tmp),]
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(tmp[tmp$time == tmp2[t,"time"] & tmp$Immu_anno == tmp2[t,"Immu_anno"],])/nrow(tmp[tmp$time == tmp2[t,"time"],])
}

position <- read.csv('data/fig1d_alluvial_plot.csv')

tmp2 <- merge(tmp2,position,by = 'time', all =T)


p1 <- ggplot(tmp2,aes(x = position,stratum = Immu_anno,alluvium = Immu_anno,
                      y = rate,fill = Immu_anno,label = Immu_anno)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(ImmuPalette),breaks = names(ImmuPalette)) +
  geom_flow(alpha = 0.3) +
  geom_stratum(width = 0.4,colour = "white", alpha = 0.7) +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        axis.text = element_text(size = 12),
        legend.position = "none")+ 
  scale_x_continuous(breaks = position$position,labels = position$time)
p1


pdf('figureS2/figS2h_alluvial_Immu.pdf', width = 10, height = 6)
p1
dev.off()

#figureS2i endo-----------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

endo <- readRDS('data/endo.RDS')

endoPalette =  c('#ef7c8e', '#fae8e0', '#b6e2d3', '#d8a7b1' ) 
names(endoPalette) <- c('AE','CE','Mural.Cell','VE')

tmp.meta = endo@meta.data

tmp = tmp.meta[,c("time","endo_anno")]
tmp$new = paste(tmp$time,tmp$endo_anno,sep = "_")
rownames(tmp) = NULL
tmp2 = tmp[!duplicated(tmp),]
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(tmp[tmp$time == tmp2[t,"time"] & tmp$endo_anno == tmp2[t,"endo_anno"],])/nrow(tmp[tmp$time == tmp2[t,"time"],])
}

position <- read.csv('data/fig1d_alluvial_plot.csv')

tmp2 <- merge(tmp2,position,by = 'time', all =T)


p1 <- ggplot(tmp2,aes(x = position,stratum = endo_anno,alluvium = endo_anno,
                      y = rate,fill = endo_anno,label = endo_anno)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(endoPalette),breaks = names(endoPalette)) +
  geom_flow(alpha = 0.3) +
  geom_stratum(width = 0.4,colour = "white", alpha = 0.7) +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        axis.text = element_text(size = 12),
        legend.position = "none")+ 
  scale_x_continuous(breaks = position$position,labels = position$time)
p1


pdf('figureS2/figS2i_alluvial_endo.pdf', width = 10, height = 6)
p1
dev.off()

