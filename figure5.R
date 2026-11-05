library(ggplot2)
library(Seurat)

load('data/mus_scRNA.Rdata')
pal11 <- c("#ffd6a5","#74B346",'#4c956c','lightgrey','lightgrey','lightgrey','lightgrey','lightgrey','lightgrey'
)

names(pal11) <-  c("Mesen","SMC",'VSMC','ESC','Epi','prolif','Endo','Immu','Eryth')

#figure5a tsne-----------------------------------------------------------------------
myeloid_fdg_coor <- scRNA@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- scRNA@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=mus_anno),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(pal11),breaks = names(pal11))+
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

ggsave(p1,filename = "figure5/fig5a_scRNA_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure5b markers-----------------------------------------------------------------------
library(scCustomize)
library(reshape2)
library(patchwork)
library(dplyr)

level = c("Endo","Epi","Immu","prolif" ,'ESC',"Mesen","SMC",'VSMC')
scRNA@meta.data[["mus_anno"]] <- factor(scRNA@meta.data[["mus_anno"]],levels = level)
Idents(scRNA) <- 'mus_anno'
markers <- FindAllMarkers(scRNA,min.pct = 0.3,logfc.threshold = 1,only.pos = T)
markers <- markers[markers$p_val_adj < 0.05,]
write.csv(markers,file = 'Table_mouse_DEGs_for_all_celltypes.csv',row.names = F)
use.gene <- markers %>% group_by(cluster) %>% do(head(., n = 4))
plot.gene <- use.gene$gene
plot.gene[17:20] <- c('Dpt', 'Igfbp6', 'Aebp1', 'Pdgfra')
plot.gene[21:24] <- c('Gas5', 'Dlk1', 'Mfap4', 'Tgfbi')
plot.gene[27:30] <- c('Des', 'Pcp4', 'Acta2', 'Mustn1')
plot.gene <- c(plot.gene,'Rgs5','Notch3')

gene_exp <- AverageExpression(scRNA,features = plot.gene,layer = 'data' )
gene_exp <- gene_exp[["RNA"]]
gene_exp <- as.data.frame(t(scale(t(gene_exp)))) 
gene_exp2 <- Percent_Expressing(scRNA,features = unique(plot.gene ))

gene_exp$gene <- row.names(gene_exp)
gene_exp2$gene <- row.names(gene_exp2)
gene_exp <- melt(gene_exp,variable.name = 'Type')
gene_exp2 <- melt(gene_exp2,variable.name = 'Type')
gene_exp3 <- merge(gene_exp,gene_exp2,by = c('gene','Type'))
gene_exp3$gene <- factor(gene_exp3$gene,levels = unique(plot.gene))


p1 = ggplot(gene_exp3,aes(x = gene , y= Type))+
  geom_point(aes(size = value.y,color = value.x ))+
  scale_color_gradient2(low='#5E3C99', high='#E66101', mid = 'white') +
  xlab("")+ylab("")+
  theme_classic()+
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))+
  labs(x = 'Gene',y = '', 
       color= 'Average Expression', size = 'Percent Expressed') 
print(p1)


pdf('figure5/fig5b_dot_plot.pdf', width = 10, height = 6)
p1 
dev.off()


#figure5c alluvial_plot-----------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

level = c("Endo","Epi","Immu","prolif" ,'ESC',"Mesen","SMC",'VSMC')
scRNA@meta.data[["mus_anno"]] <- factor(scRNA@meta.data[["mus_anno"]],levels = level)
Idents(scRNA) <- 'mus_anno'

tmp.meta = scRNA@meta.data

tmp = tmp.meta[,c("time","mus_anno")]
tmp$new = paste(tmp$time,tmp$mus_anno,sep = "_")
rownames(tmp) = NULL
tmp2 = tmp[!duplicated(tmp),]
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(tmp[tmp$time == tmp2[t,"time"] & tmp$mus_anno == tmp2[t,"mus_anno"],])/nrow(tmp[tmp$time == tmp2[t,"time"],])
}

position <- read.csv('data/mouse_alluvial_plot.csv')

tmp2 <- merge(tmp2,position,by = 'time', all =T)


p1 <- ggplot(tmp2,aes(x = position,stratum = mus_anno,alluvium = mus_anno,
                      y = rate,fill = mus_anno,label = mus_anno)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(pal11),breaks = names(pal11)) +
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



pdf('figure5/fig5c_alluvial_plot.pdf', width = 10, height = 6)
p1
dev.off()



#figure5d fa2-----------------------------------------------------------------------
mus_fa2 <-  readRDS('data/mus_fa2.RDS')
mus_fa2 <- subset(mus_fa2,cells = row.names(mus_fa2@reductions[["fa2"]]@cell.embeddings))

myeloid_fdg_coor <- mus_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- mus_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=mus_anno),shape = 21,stroke = 0.2,size = 3)+
  theme_bw()+
  scale_fill_manual(values = as.character(pal11),breaks = names(pal11))+
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
ggsave(p1,filename = "figure5/fig5d_mus_fa2.pdf",width = 6, height = 6, units = "in", device='pdf')


#figure3e cytotrace---------------------------------------------------------------------
mus_fa2 <-  readRDS('data/mus_fa2.RDS')
mus_fa2 <- subset(mus_fa2,cells = row.names(mus_fa2@reductions[["fa2"]]@cell.embeddings))

load('data/mus_CytoTRACE.Rdata')
cytodf <- as.data.frame(results[["CytoTRACE"]])
mus_cyto <- AddMetaData(mus_fa2,cytodf,col.name = "CytoTrace")

myeloid_fdg_coor <- mus_cyto@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- mus_cyto@meta.data
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
  )
p1
ggsave(p1,filename = "figure5/fig5e_cyto_legend.pdf",width = 6, height = 6, units = "in", device='pdf')

p1<-ggplot(myeloid)+
  geom_boxplot(aes(x = reorder(mus_anno, -CytoTrace), y = CytoTrace,fill=mus_anno))+
  scale_fill_manual(values = as.character(pal11),breaks = names(pal11))+
  theme_classic()+
  NoLegend()+
  theme(axis.title.x=element_blank())

p1
ggsave(p1,filename = "figure5/fig5e_cyto_box.pdf",width = 6, height = 6, units = "in", device='pdf')


#figure5f MP-----------------------------------------------------------------------
mus_fa2 <- readRDS(file = 'data/mus_fa2_MP.RDS')

myeloid_fdg_coor <- mus_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- mus_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)
colnames(myeloid)[19:36] <- paste0('MP',colnames(myeloid)[19:36])

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=MP2),shape = 21,stroke = 0.2,size = 3)+
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
ggsave(p1,filename = "figure5/fig5f_mus_MP2.pdf",width = 6, height = 6, units = "in", device='pdf')

p2<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=MP12),shape = 21,stroke = 0.2,size = 3)+
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

p2
ggsave(p2,filename = "figure5/fig5f_mus_MP12.pdf",width = 6, height = 6, units = "in", device='pdf')

p3<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=MP9),shape = 21,stroke = 0.2,size = 3)+
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

p3
ggsave(p3,filename = "figure5/fig5f_mus_MP9.pdf",width = 6, height = 6, units = "in", device='pdf')

p4<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=MP10),shape = 21,stroke = 0.2,size = 3)+
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

p4
ggsave(p4,filename = "figure5/fig5f_mus_MP10.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure5h dorothea-----------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(tibble)
library(pheatmap)

load('data/dorothea_mouse.Rdata')
mesen_fa2 <- ScaleData(mesen_fa2)

viper_scores_df <- GetAssayData(mesen_fa2, slot = "scale.data", assay = "dorothea") %>%
  data.frame(check.names = F) %>% t()
viper_scores_df[1:4,1:4]

CellsClusters <- data.frame(cell = names(Idents(mesen_fa2)), 
                            cell_type = as.character(Idents(mesen_fa2)),
                            check.names = F)
head(CellsClusters)

viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)


summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

head(summarized_viper_scores)

highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(147, var) %>%
  distinct(tf)
highly_variable_tfs <- rbind(highly_variable_tfs,'Twist1')

summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

summarized_viper_scores_df[1:3,1:4]
summarized_viper_scores_df <- t(summarized_viper_scores_df)
summarized_viper_scores_df <- summarized_viper_scores_df[,c('VSMC','Mesen','SMC')]

palette_length = 100
my_color = colorRampPalette(c("#3e5aa8", "white","#ff6361"))(100)

p1 <- pheatmap(summarized_viper_scores_df,fontsize=14, 
               fontsize_row = 10, cluster_cols = F,
               color=my_color, breaks = seq(-0.5, 0.5, length.out = 101), 
               main = "DoRothEA (ABC)", angle_col = 45,
               treeheight_col = 0,  border_color = NA) 


pdf('figure5/fig5h_dorothea.pdf', width = 8, height = 8)
p1
dev.off()
