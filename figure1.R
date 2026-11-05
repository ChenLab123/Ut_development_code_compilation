library(ggplot2)
library(Seurat)

load('data/JoinLayers_scRNA.rdata')

pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff"
)

names(pal11) <-  c("Endo","Epi.MD", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen")

group_colors<-c( '#ff6361', '#bc5090', '#58508d' )
names(group_colors) <- c('Prenatal stage','Child-bearing stage','Postmenopause stage')


#figure1b-----------------------------------------------------------------
library(plot1cell)

colnames(combined@meta.data)  

combined@reductions[["umap"]] <- combined@reductions[["umap.harmony"]]
Idents(combined) <- 'harmony_anno27'

circ_data <- prepare_circlize_data(combined, scale = 0.8 )

pdf('figure1/fig1b_circlize_plot.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = pal11 ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
add_track(circ_data, group = "Stage", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
dev.off()

#figure1c-----------------------------------------------------------------
library(scCustomize)
library(reshape2)
library(patchwork)
library(dplyr)

level = c("Endo","Epi.MD","Eryth","ESC","Immu","Mesen","ESR1.Mesen","THY1.Mesen",
          "NCC" ,"prolif" ,"SMC")
combined@meta.data[["harmony_anno27"]] <- factor(combined@meta.data[["harmony_anno27"]],levels = level)
Idents(combined) <- 'harmony_anno27'
markers <- FindAllMarkers(combined,min.pct = 0.3,logfc.threshold = 1,only.pos = T)
markers <- markers[markers$p_val_adj < 0.05,]
write.csv(markers,file = 'TableS2_DEGs_for_all_celltypes.csv',row.names = F)
use.gene <- markers %>% group_by(cluster) %>% do(head(., n = 2))

plot.gene <- use.gene$gene
plot.gene <- append(plot.gene, "ESR1", after = which(plot.gene == "COL14A1"))
plot.gene <- append(plot.gene, "THY1", after = which(plot.gene == "ELN"))

gene_exp <- AverageExpression(combined,features = plot.gene,layer = 'data' )
gene_exp <- gene_exp[["RNA"]]
gene_exp <- as.data.frame(t(scale(t(gene_exp)))) 
gene_exp2 <- Percent_Expressing(combined,features = plot.gene )

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

cluster <- c("Endo","Epi.MD","Eryth","ESC","Immu","Mesen","ESR1.Mesen","THY1.Mesen",
            "NCC" ,"prolif" ,"SMC")
cluster <- factor(cluster,levels = rev(level)) 
type <- 'Bar'
value <- 1
data <- data.frame(cluster,type,value)

p2 <- ggplot(data, aes(fill=cluster, y=value, x=type)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = pal11,breaks = names(pal11))+
  theme_void()+
  theme(legend.position = 'none')


pdf('figure1/fig1c_dot_plot.pdf', width = 10, height = 6)
p2 + p1 + plot_layout(widths = c(1, 20))
dev.off()

#figure1d-----------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

tmp.meta = combined@meta.data

tmp = tmp.meta[,c("time","harmony_anno27")]
tmp$new = paste(tmp$time,tmp$harmony_anno27,sep = "_")
rownames(tmp) = NULL
tmp2 = tmp[!duplicated(tmp),]
tmp2$rate = 0
for(t in rownames(tmp2)){
    tmp2[t,"rate"] = nrow(tmp[tmp$time == tmp2[t,"time"] & tmp$harmony_anno27 == tmp2[t,"harmony_anno27"],])/nrow(tmp[tmp$time == tmp2[t,"time"],])
  }

position <- read.csv('data/fig1d_alluvial_plot.csv')

tmp2 <- merge(tmp2,position,by = 'time', all =T)


p1 <- ggplot(tmp2,aes(x = position,stratum = harmony_anno27,alluvium = harmony_anno27,
                y = rate,fill = harmony_anno27,label = harmony_anno27)) +
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


pdf('figure1/fig1d_alluvial_plot.pdf', width = 10, height = 6)
p1
dev.off()

