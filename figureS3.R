library(ggplot2)
library(Seurat)

#figureS3b tsne-----------------------------------------------------------------
mesen <- readRDS('data/mesen.RDS')
mesen <- subset(mesen,subset = mesen_anno %in% c('Mesen3','Mesen7'),invert = T)
mesen <- subset(mesen,cells =  c('Ut26_GTTCCGTTCTTCGCTG-1','GTTCCGTGTAACTTCG-1_2_1_1'),invert = T)
mesen$subanno <- ifelse(mesen$mesen_anno %in% c('Mesen0','Mesen2','Mesen6'),mesen$mesen_anno,'Others')

mesenPalette =  c("#a0c4ff", "#ffd6a5", "#00c2ff") 
names(mesenPalette) <- c("ESR1.Mesen", "Mesen","THY1.Mesen")

myeloid_fdg_coor <- mesen@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- mesen@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=harmony_anno27),shape = 21,stroke = 0.2,size = 3)+
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
ggsave(p1,filename = "figureS3/figS3b_mesen_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS3c umap-----------------------------------------------------------------
mesen <- readRDS('data/mesen.RDS')
mesen$subanno <- ifelse(mesen$mesen_anno %in% c('Mesen0','Mesen2','Mesen6'),mesen$mesen_anno,'Others')

mesendf <- mesen$subanno

mesenPalette =  c("#ffadad",'#2c6e49', '#4c956c','#fefee3', '#d68c45', '#ffc9b9','#ac6f82') 
names(mesenPalette) <- c("ESC","DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6','Others')


SEM <- readRDS('data/SEM_subanno.RDS')
SEM$subanno <- ifelse(SEM$subanno %in% c('DES.SMC','VSMC'),SEM$subanno,as.character(SEM$harmony_anno27))
SEM <- subset(SEM ,subset = subanno == 'SMC',invert = T)
SEMdf <- SEM$subanno

load('data/JoinLayers_scRNA.rdata')

combined <- AddMetaData(combined,mesendf,col.name = 'mesen_anno')
combined <- AddMetaData(combined,SEMdf,col.name = 'SEM_anno')
combined$subanno <- ifelse(is.na(combined$mesen_anno),combined$SEM_anno,combined$mesen_anno)

myeloid_fdg_coor <- combined@reductions[["umap.harmony"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- combined@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)
myeloid <- myeloid[!is.na(myeloid$subanno),]

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
ggsave(p1,filename = "figureS_260415/figS_SEMsub_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS3f -----------------------------------------------------------------
load('data/ut_12_008.Rdata')
qc12df <- object@meta.data
qc12df$time <- 'PCW8'
load('data/ut_15_008.Rdata')
qc15df <- object@meta.data
qc15df$time <- 'PCW11'

qcdf <- rbind(qc12df,qc15df)

qcdf <- qcdf[,c("nCount_Spatial.008um","nFeature_Spatial.008um","time")]
qcdf$time <- factor(qcdf$time,levels = c('PCW8','PCW11'))

sample_size = qcdf %>% group_by(time) %>% summarize(num=n())
sample_size$myaxis = paste0(sample_size$time, "\n", "n=", sample_size$num)

p1 <- ggplot(qcdf, aes(x=time, y=nCount_Spatial.008um, fill=time)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, color="grey", alpha=0.2,outlier.shape = NA) +
  scale_fill_manual(values = c("#440154FF","#3E4A89FF")) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("") +
  xlab("") +
  coord_cartesian(ylim = c(0, 3000))+ 
  scale_x_discrete(labels = sample_size$myaxis)

p2 <- ggplot(qcdf, aes(x=time, y=nFeature_Spatial.008um, fill=time)) +
  geom_violin(scale = "width") +
  geom_boxplot(width=0.1, color="grey", alpha=0.2,outlier.shape = NA) +
  scale_fill_manual(values = c("#440154FF","#3E4A89FF")) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("") +
  xlab("") +
  coord_cartesian(ylim = c(0, 2500))+ 
  scale_x_discrete(labels = sample_size$myaxis)

p1/p2

ggsave(p1/p2,filename = "figureS3/figS3f_sp_QC.pdf",width = 3, height = 6, units = "in", device='pdf')

#figureS3g -----------------------------------------------------------------
pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff",'#ffc9b9'
)

names(pal11) <-  c("Endo","Epi.MD", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen",'Mesen6')

myRCTD <- readRDS(file = 'data/scRNA_27_12myRCTD_subanno6_2.rds')
load('data/ut_12_008.Rdata')

rctd <- myRCTD@results$results_df

object<- subset(object,cells = row.names(rctd))
object <- AddMetaData(object,rctd)
object <- RunTSNE(object, reduction = "pca.008um", reduction.name = "tsne.008um", dims = 1:30)

myeloid_fdg_coor <- object@reductions[["tsne.008um"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- object@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=first_type),shape = 21,stroke = 0.2,size = 3)+
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
  )

p1
ggsave(p1,filename = "figureS3/figS3g_PCW8_withlegend_2.pdf",width = 6, height = 6, units = "in", device='pdf')


myRCTD <- readRDS(file = 'data/scRNA_27_15myRCTD_subanno6_2.rds')
load('data/ut_15_008.Rdata')

rctd <- myRCTD@results$results_df

object<- subset(object,cells = row.names(rctd))
object <- AddMetaData(object,rctd)
object <- RunTSNE(object, reduction = "pca.008um", reduction.name = "tsne.008um", dims = 1:30)

myeloid_fdg_coor <- object@reductions[["tsne.008um"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- object@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=first_type),shape = 21,stroke = 0.2,size = 3)+
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
  )

p1
ggsave(p1,filename = "figureS3/figS3g_PCW11_withlegend_2.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS3h dim plot-----------------------------------------------------------------------
pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff",'#ffc9b9'
)

names(pal11) <-  c("Endo","Epi.MD", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen",'Mesen6')

myRCTD <- readRDS(file = 'data/scRNA_27_12myRCTD_subanno6_2.rds')
load('data/ut_12_008.Rdata')

rctd <- myRCTD@results$results_df

coords <- GetTissueCoordinates(object)[, 1:2]

spdf <- merge(coords,rctd,by = 0)


p2<-ggplot(spdf)+
  geom_point(aes(x = x, y = y,fill=first_type),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
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
  )
print(p2)
ggsave(p2,filename = "figureS3/FigS3h_PCW8_tolta_2.pdf",width = 8, height = 5, units = "in", device='pdf')


myRCTD <- readRDS(file = 'data/scRNA_27_15myRCTD_subanno6_2.rds')
load('data/ut_15_008.Rdata')

rctd <- myRCTD@results$results_df

coords <- GetTissueCoordinates(object)[, 1:2]

spdf <- merge(coords,rctd,by = 0)


p2<-ggplot(spdf)+
  geom_point(aes(x = x, y = y,fill=first_type),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
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
  )
print(p2)
ggsave(p2,filename = "figureS3/FigS3h_PCW11_tolta_2.pdf",width = 8, height = 5, units = "in", device='pdf')

