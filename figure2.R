library(ggplot2)
library(Seurat)

#figure2c tsne-----------------------------------------------------------------
mesen <- readRDS('data/mesen.RDS')
mesen <- subset(mesen,subset = mesen_anno %in% c('Mesen3','Mesen7'),invert = T)
mesen <- subset(mesen,cells =  c('Ut26_GTTCCGTTCTTCGCTG-1','GTTCCGTGTAACTTCG-1_2_1_1'),invert = T)
mesen$subanno <- ifelse(mesen$mesen_anno %in% c('Mesen0','Mesen2','Mesen6'),mesen$mesen_anno,'Others')

mesenPalette =  c( '#fefee3', '#d68c45', '#ffc9b9','#ac6f82') 
names(mesenPalette) <- c('Mesen0','Mesen2','Mesen6','Others')

myeloid_fdg_coor <- mesen@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- mesen@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=subanno),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(mesenPalette),breaks = names(mesenPalette))+
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
ggsave(p1,filename = "figure2/fig2c_mesen_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figure2c vln-----------------------------------------------------------------
library(scCustomize)
library(reshape2)
library(patchwork)
library(dplyr)

mesenPalette =  c( '#fefee3', '#d68c45', '#ffc9b9','#ac6f82') 
names(mesenPalette) <- c('Mesen0','Mesen2','Mesen6','Others')

mesen <- readRDS('data/mesen.RDS')
mesen$subanno <- ifelse(mesen$mesen_anno %in% c('Mesen0','Mesen2','Mesen6'),mesen$mesen_anno,'Others')
mesen@meta.data[["subanno"]] <- factor(mesen@meta.data[["subanno"]],levels = c('Mesen6','Mesen0','Others','Mesen2'))
Idents(mesen) <- 'subanno'

plot.gene <- c("ACTA2","TAGLN","ALDH1A1","ENPP2","DPP6","COL25A1",'DCN','COL1A2',
               'OGN','SFRP2','SULT1E1','ESR1', 'ELN', 'THY1', 'SFRP1', 'MYH11')
plot.gene <- plot.gene[order(plot.gene)]

p1 <- VlnPlot(mesen,plot.gene,stack = T,fill.by = 'ident')+
  scale_fill_manual(values = as.character(mesenPalette),breaks = names(mesenPalette))

cluster <- c('Mesen0','Mesen2','Mesen6','Others')
cluster <- factor(cluster,levels = c('Mesen2','Others','Mesen0','Mesen6')) 
type <- 'Bar'
value <- 1
data <- data.frame(cluster,type,value)

p2 <- ggplot(data, aes(fill=cluster, y=value, x=type)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = mesenPalette,breaks = names(mesenPalette))+
  theme_void()+
  theme(legend.position = 'none')

pdf('figure2/fig2c_mesen_vln.pdf', width = 10, height = 4)
p2 + p1 + plot_layout(widths = c(1, 20))
dev.off()

#figure2d alluvial-----------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

mesen <- readRDS('data/mesen.RDS')
mesen$subanno <- ifelse(mesen$mesen_anno %in% c('Mesen0','Mesen2','Mesen6'),mesen$mesen_anno,'Others')
mesen <- subset(mesen,subset = time %in% c('Adult33', 'Adult55', 'Adult65', 'Adult71'),invert = T)

tmp.meta = mesen@meta.data

tmp = tmp.meta[,c("time","subanno")]
tmp$new = paste(tmp$time,tmp$subanno,sep = "_")
rownames(tmp) = NULL
tmp2 = tmp[!duplicated(tmp),]
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(tmp[tmp$time == tmp2[t,"time"] & tmp$subanno == tmp2[t,"subanno"],])/nrow(tmp[tmp$time == tmp2[t,"time"],])
}

position <- read.csv('data/fig1d_alluvial_plot.csv')

tmp2 <- merge(tmp2,position,by = 'time', all.x =T)


p1 <- ggplot(tmp2,aes(x = position,stratum = subanno,alluvium = subanno,
                      y = rate,fill = subanno,label = subanno)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(mesenPalette),breaks = names(mesenPalette)) +
  geom_flow(alpha = 0.3) +
  geom_stratum(width = 0.4,colour = "#5E3C99", alpha = 0.7) +
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

pdf('figure2/fig2d_mesen_alluvial.pdf', width = 5, height = 6)
p1
dev.off()

#figure2e-----------------------------------------------------------------
library(ggrepel)
SEM <- readRDS('data/SEM_subanno.RDS')
gene_exp <- AverageExpression(SEM,layer = 'data' ,group.by = 'subanno')
gene_exp <- gene_exp[["RNA"]]
gene_exp <- as.data.frame(gene_exp)
gene_exp <- log2(gene_exp+1)

gene_exp$Mesen6vsDES.SMC <- ifelse((gene_exp$Mesen6- gene_exp$DES.SMC > 0.5),"Up",
                                   ifelse((gene_exp$Mesen6 - gene_exp$DES.SMC < -0.5),"Down","Others"))

genelabel <- c('ACTA2', 'MYH11', 'TAGLN', 'DES', 'RGS5', 'NOTCH3', 'PDGFRB', 'IGFBP3', 'CXCL12', 'NCAM1')
gene_exp$label = ifelse(row.names(gene_exp) %in% genelabel, row.names(gene_exp), "")

p1 <- ggplot(gene_exp, aes(x = DES.SMC, y = Mesen6,fill = Mesen6vsDES.SMC)) +
  geom_point(shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0)) +
  geom_abline(
    intercept = 0,slope = 1,linetype = "dashed",color = "black",linewidth = 0.2) +
  geom_abline(intercept = 0.5,slope = 1,linetype = "dotted",color = "black",linewidth = 0.2) +
  geom_abline(intercept = -0.5,slope = 1,linetype = "dotted",color = "black",linewidth = 0.2) +
  coord_fixed(ratio = 1) +
  labs(x = "DES.SMC", y = "Mesen6", title = "Mesen6 vs DES.SMC") +
  theme_bw()+
  scale_fill_manual(values = c("Others" = "grey","Up" = "#ffc9b9","Down" ="#2c6e49"))+
  geom_text_repel(aes(label = label), size = 3, box.padding = unit(0.5, "lines"), 
                  max.overlaps = Inf,min.segment.length = 0, force = 20,force_pull = 0.5, 
                  point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)+
  NoLegend()
p1
ggsave(p1,filename = "figure2/Fig2e_DES.SMC.pdf",width = 8, height = 8, units = "in", device='pdf')


gene_exp$Mesen6vsVSMC <- ifelse((gene_exp$Mesen6 - gene_exp$VSMC > 0.5),"Up",
                                ifelse((gene_exp$Mesen6 - gene_exp$VSMC < -0.5),"Down","Others"))

genelabel <- c('ACTA2', 'MYH11', 'TAGLN', 'DES', 'RGS5', 'NOTCH3', 'PDGFRB', 'IGFBP3', 'CXCL12', 'NCAM1')
gene_exp$label = ifelse(row.names(gene_exp) %in% genelabel, row.names(gene_exp), "")

p1 <-ggplot(gene_exp, aes(x = VSMC, y = Mesen6,fill = Mesen6vsVSMC)) +
  geom_point(shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0)) +
  geom_abline(
    intercept = 0,slope = 1,linetype = "dashed",color = "black",linewidth = 0.2) +
  geom_abline(intercept = 0.5,slope = 1,linetype = "dotted",color = "black",linewidth = 0.2) +
  geom_abline(intercept = -0.5,slope = 1,linetype = "dotted",color = "black",linewidth = 0.2) +
  coord_fixed(ratio = 1) +
  labs(x = "VSMC", y = "Mesen6", title = "Mesen6 vs VSMC") +
  theme_bw()+
  scale_fill_manual(values = c("Others" = "grey","Up" = "#ffc9b9","Down" ="#4c956c"))+
  geom_text_repel(aes(label = label), size = 3, box.padding = unit(0.5, "lines"), 
                  max.overlaps = Inf,min.segment.length = 0, force = 20,force_pull = 0.5, 
                  point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)+
  NoLegend()
p1
ggsave(p1,filename = "figure2/Fig2e_VSMC.pdf",width = 8, height = 8, units = "in", device='pdf')



#figure2g-----------------------------------------------------------------
library(ggsignif)
pal7 <- c( '#fefee3', '#d68c45', '#ffc9b9','#ac6f82',"#ffa7dc","#debcff","#caffbf" )
names(pal7) <-  c('Mesen0','Mesen2','Mesen6','Others','prolif','Endo','Epi.MD')

load('data/ut_12_008.Rdata')
genelist <- c('OGN', 'SULT1E1')

celldf <- object@assays$Spatial.008um$data[genelist,]
celldf <- as.data.frame(t(celldf))

coords <- GetTissueCoordinates(object)[, 1:2]

spdf <- merge(coords,celldf,by = 0)
spdf <- spdf[spdf$x > 4000 & spdf$x < 11000 & spdf$y < 16000 & spdf$y > 9000,]

col <- c("#1A0000", "#8B0000", "#FF2200", "#FF5500", "#FFAA00", "#FFDD00")


pdf('figure2/fig2g_PCW8_markers.pdf',width = 4, height = 4)
for (i in genelist) {
  
  p1 <- ggplot(spdf)+
    geom_point(aes_string(x = 'x', y = 'y',fill=i),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank()
    )+
    scale_fill_gradientn(colors = col) 
  print(p1)
}
dev.off()


myRCTD <- readRDS(file = 'data/scRNA_27_12myRCTD_260424.rds')

rctd <- myRCTD@results$results_df

coords <- GetTissueCoordinates(object)[, 1:2]

spdf <- merge(coords,rctd,by = 0)
spdf <- spdf[spdf$x > 4000 & spdf$x < 11000 & spdf$y < 16000 & spdf$y > 9000,]

p1<-ggplot(spdf)+
  geom_point(aes(x = x, y = y,fill=first_type),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  scale_fill_manual(values = as.character(pal7),breaks = names(pal7))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )
print(p1)

ggsave(p1,filename = "figure2/Fig2g_PCW8.pdf",width = 4, height = 5, units = "in", device='pdf')

out1 <- spdf[spdf$x > 9500 & spdf$x < 10000 & spdf$y > 15400 & spdf$y < 16000,]
out1$type <- 'out1'
out2 <- spdf[spdf$x > 7000 & spdf$x < 7500 & spdf$y > 15000 & spdf$y < 15500,]
out2$type <- 'out2'
out3 <- spdf[spdf$x > 10500 & spdf$x < 11000 & spdf$y > 11000 & spdf$y < 11600,]
out3$type <- 'out3'
in1 <- spdf[spdf$x > 6500 & spdf$x < 7000 & spdf$y > 13500 & spdf$y < 14000,]
in1$type <- 'in1'
in2 <- spdf[spdf$x > 8000 & spdf$x < 8500 & spdf$y > 11500 & spdf$y < 12000,]
in2$type <- 'in2'
in3 <- spdf[spdf$x > 8500 & spdf$x < 9000 & spdf$y > 13500 & spdf$y < 14000,]
in3$type <- 'in3'

typedf <- rbind(out1,out2,out3,in1,in2,in3)
typedf <- typedf[,c("Row.names","type")]

spdf <- merge(spdf,typedf,by = 'Row.names',all = T)

typepal <- c(rep('#00c2ff',3),rep('#a0c4ff',3))
names(typepal) <- c("out1","out2", "out3","in1","in2","in3")

p1<-ggplot(spdf)+
  geom_point(aes(x = x, y = y,fill=type),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  scale_fill_manual(values = as.character(typepal),breaks = names(typepal))+
  theme_bw()+
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

ggsave(p1,filename = "figure2/Fig2g_PCW8_split.pdf",width = 4, height = 5, units = "in", device='pdf')

genedf <- object@assays$Spatial.008um$data
genedf <- t(as.matrix(genedf)[genelist,typedf$Row.names])
genedf <- merge(genedf,typedf,by.x = 0,by.y = 'Row.names')
genedf <- genedf[,-1]

group_mean <- aggregate(.~ type, data = genedf, mean) 
group_mean$Type <- c('In','In','In','Out','Out','Out')
group_mean$Type <- factor(group_mean$Type,levels = c('Out','In'))

p1 <- ggplot(group_mean,aes(Type,OGN))+
  geom_bar(aes(fill=Type),stat="summary",fun=mean,position="dodge")+
  geom_jitter(color="black",size=1,width=0.2,alpha=0.8)+
  scale_fill_manual(values=c("#00c2ff","#a0c4ff"))+
  theme_classic()+
  geom_signif(comparisons = list(c('Out','In')),map_signif_level = T,test = 't.test')

p2 <- ggplot(group_mean,aes(Type,SULT1E1))+
  geom_bar(aes(fill=Type),stat="summary",fun=mean,position="dodge")+
  geom_jitter(color="black",size=1,width=0.2,alpha=0.8)+
  scale_fill_manual(values=c("#00c2ff","#a0c4ff"))+
  theme_classic()+
  geom_signif(comparisons = list(c('Out','In')),map_signif_level = T,test = 't.test')


p1+p2

ggsave((p1+p2),filename = "figure2/Fig2g_PCW8_genebar.pdf",width = 8, height = 8, units = "in", device='pdf')



load('data/ut_15_008.Rdata')
genelist <- c('OGN', 'SULT1E1')

celldf <- object@assays$Spatial.008um$data[genelist,]
celldf <- as.data.frame(t(celldf))

coords <- GetTissueCoordinates(object)[, 1:2]

spdf <- merge(coords,celldf,by = 0)
spdf <- spdf[spdf$x > -30000 & spdf$x < -15000 & spdf$y < 24000 & spdf$y > 16000,]

col <- c("#1A0000", "#8B0000", "#FF2200", "#FF5500", "#FFAA00", "#FFDD00")

pdf('figure2/fig2g_PCW11_markers.pdf',width = 4, height = 4)
for (i in genelist) {
  
  p1 <- ggplot(spdf)+
    geom_point(aes_string(x = 'x', y = 'y',fill=i),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank()
    )+
    scale_fill_gradientn(colors = col) 
  print(p1)
}
dev.off()

myRCTD <- readRDS(file = 'data/scRNA_27_15myRCTD_260424.rds')

rctd <- myRCTD@results$results_df

coords <- GetTissueCoordinates(object)[, 1:2]

spdf <- merge(coords,rctd,by = 0)
spdf <- spdf[spdf$x > -30000 & spdf$x < -15000 & spdf$y < 24000 & spdf$y > 16000,]


p1<-ggplot(spdf)+
  geom_point(aes(x = x, y = y,fill=first_type),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  scale_fill_manual(values = as.character(pal7),breaks = names(pal7))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank()
  )
print(p1)

ggsave(p1,filename = "figure2/Fig2g_PCW11.pdf",width = 5, height = 5, units = "in", device='pdf')


out1 <- spdf[spdf$x > -27500 & spdf$x < -27000 & spdf$y > 22500 & spdf$y < 23000,]
out1$type <- 'out1'
out2 <- spdf[spdf$x > -29500 & spdf$x < -29000 & spdf$y > 21000 & spdf$y < 21500,]
out2$type <- 'out2'
out3 <- spdf[spdf$x > -29500 & spdf$x < -29000 & spdf$y > 16800 & spdf$y < 17300,]
out3$type <- 'out3'
in1 <- spdf[spdf$x > -27500 & spdf$x < -27000 & spdf$y > 21000 & spdf$y < 21500,]
in1$type <- 'in1'
in2 <- spdf[spdf$x > -29500 & spdf$x < -29000 & spdf$y > 20000 & spdf$y < 20500,]
in2$type <- 'in2'
in3 <- spdf[spdf$x > -28500 & spdf$x < -28000 & spdf$y > 18000 & spdf$y < 18500,]
in3$type <- 'in3'

typedf <- rbind(out1,out2,out3,in1,in2,in3)
typedf <- typedf[,c("Row.names","type")]

spdf <- merge(spdf,typedf,by = 'Row.names',all = T)

typepal <- c(rep('#00c2ff',3),rep('#a0c4ff',3))
names(typepal) <- c("out1","out2", "out3","in1","in2","in3")

p1<-ggplot(spdf)+
  geom_point(aes(x = x, y = y,fill=type),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  scale_fill_manual(values = as.character(typepal),breaks = names(typepal))+
  theme_bw()+
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

ggsave(p1,filename = "figure2/Fig2g_PCW11_split.pdf",width = 5, height = 5, units = "in", device='pdf')

genedf <- object@assays$Spatial.008um$data
genedf <- t(as.matrix(genedf)[genelist,typedf$Row.names])
genedf <- merge(genedf,typedf,by.x = 0,by.y = 'Row.names')
genedf <- genedf[,-1]

group_mean <- aggregate(.~ type, data = genedf, mean) 
group_mean$Type <- c('In','In','In','Out','Out','Out')
group_mean$Type <- factor(group_mean$Type,levels = c('Out','In'))

p1 <- ggplot(group_mean,aes(Type,OGN))+
  geom_bar(aes(fill=Type),stat="summary",fun=mean,position="dodge")+
  geom_jitter(color="black",size=1,width=0.2,alpha=0.8)+
  scale_fill_manual(values=c("#00c2ff","#a0c4ff"))+
  theme_classic()+
  geom_signif(comparisons = list(c('Out','In')),map_signif_level = T,test = 't.test')
p2 <- ggplot(group_mean,aes(Type,SULT1E1))+
  geom_bar(aes(fill=Type),stat="summary",fun=mean,position="dodge")+
  geom_jitter(color="black",size=1,width=0.2,alpha=0.8)+
  scale_fill_manual(values=c("#00c2ff","#a0c4ff"))+
  theme_classic()+
  geom_signif(comparisons = list(c('Out','In')),map_signif_level = T,test = 't.test')

p1+p2

ggsave((p1+p2),filename = "figure2/Fig2g_PCW11_genebar.pdf",width = 8, height = 8, units = "in", device='pdf')

