library(scCustomize)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)
library(Seurat)

#figureS4a mesen fa2----------------------------------------------------------------------------
mesen_fa2 <- readRDS("data/mesen_fa2.RDS")

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

myeloid_fdg_coor <- mesen_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- mesen_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)
myeloid$subanno <- ifelse(myeloid$subanno %in% c('RERGL.VSMC','STEAP4.VSMC','RGS6.VSMC'),
                          'VSMC',as.character(myeloid$subanno)) 
library(dplyr)
myeloid <- myeloid %>%
  arrange(subanno == 'Mesen6')  

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
ggsave(p1,filename = "figureS4/figS4a_mesen_fa2.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS4b VSMC fa2----------------------------------------------------------------------------
VSMC_fa2 <-  readRDS('data/VSMC_fa2.RDS')

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

myeloid_fdg_coor <- VSMC_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- VSMC_fa2@meta.data
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
ggsave(p1,filename = "figureS4/figS4b_VSMC_fa2.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS4c fa2---------------------------------------------------------------------
sp_fa2 <-  readRDS('data/sp_fa2.RDS')

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')

myeloid_fdg_coor <- sp_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- sp_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=first_type),shape = 21,stroke = 0.2,size = 3)+
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
ggsave(p1,filename = "figureS4/figS4c_sp_fa2.pdf",width = 6, height = 6, units = "in", device='pdf')

#figureS4de cytotrace---------------------------------------------------------------------
sp_fa2 <- readRDS('data/sp_fa2.RDS')
load('data/sp_CytoTRACE.Rdata')
cytodf <- as.data.frame(results[["CytoTRACE"]])
SMC_cyto <- AddMetaData(sp_fa2,cytodf,col.name = "CytoTrace")

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
  )
p1
ggsave(p1,filename = "figureS4/figS4d_cyto.pdf",width = 6, height = 6, units = "in", device='pdf')

SEMPalette =  c("#ffadad","#a0c4ff", "#ffd6a5",  "#74B346", "#00c2ff",
                '#2c6e49', '#4c956c', '#fefee3', '#d68c45', '#ffc9b9') 
names(SEMPalette) <- c("ESC", "ESR1.Mesen", "Mesen","SMC","THY1.Mesen",
                       "DES.SMC","VSMC",'Mesen0','Mesen2','Mesen6')
p1<-ggplot(myeloid)+
  geom_boxplot(aes(x = reorder(first_type, -CytoTrace), y = CytoTrace,fill=first_type))+
  scale_fill_manual(values = as.character(SEMPalette),breaks = names(SEMPalette))+
  theme_classic()+
  NoLegend()+
  theme(axis.title.x=element_blank())

p1
ggsave(p1,filename = "figureS4/figS4e_cyto_box.pdf",width = 6, height = 6, units = "in", device='pdf')
