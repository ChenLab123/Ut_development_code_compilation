#figureS6b matrix_stiffness_gene_modules----------------------------------------------------------------------------
library(Seurat)
library(UCell)
library(dplyr)
library(reshape2)
library(tidyr)
library(patchwork)


matrix_stiffness_gene_modules <- list(ECM_structural_core = c("ACAN", "APCS", "ASPN", "BCAN", "BGN", "C1QL1", "C1QL2", "C1QL3", "C1QL4", "C1QTNF1", "C1QTNF2", "C1QTNF3", "C1QTNF5", "C1QTNF6", "C1QTNF7", "C1QTNF8", "C1QTNF9", "C1QTNF9B", "CCBE1", "CHADL", "CILP", "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", "COL19A1", "COL1A2", "COL20A1", "COL21A1", "COL22A1", "COL23A1", "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1", "COL2A1", "COL3A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A6", "COL5A1", "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6", "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COLEC10", "COLEC11", "COLEC12", "COLGALT1", "COLQ", "COMP", "DCN", "DMP1", "DPT", "ECM1", "ECM2", "EFEMP1", "EFEMP2", "ELN", "EMILIN1", "EMILIN2", "EMILIN3", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FBN3", "FGA", "FGB", "FGG", "FMOD", "FN1", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6", "HAPLN1", "HAPLN2", "HAPLN3", "HAPLN4", "HMCN1", "HMCN2", "HSPG2", "IBSP", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", "LAMB4", "LAMC1", "LAMC2", "LAMC3", "LUM", "MATN1", "MATN2", "MATN3", "MATN4", "MFAP2", "MFAP4", "MFAP5", "MGP", "MMP23B", "NID1", "NID2", "NTN1", "NTN3", "NTN4", "OGN", "POSTN", "PRELP", "SFTPA", "SFTPA2", "SFTPD", "SLIT2", "SMOC1", "SMOC2", "SPARC", "SPON1", "SPP1", "SPP2", "SRPX", "SRPX2", "THBS1", "THBS2", "THBS3", "THBS4", "TINAG", "TINAGL1", "TNC", "TNN", "TNR", "TNXA", "TNXB", "VIT", "VTN", "VWA1", "VWA2", "VWF"),
                                      ECM_crosslinking_maturation = c("CRTAP", "FKBP10", "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "P3H2", "P3H3", "P4HA1", "P4HA2", "P4HA3", "P4HB", "PCOLCE", "PCOLCE2", "PLOD3", "PXDN", "SERPINH1", "TGM2", "TGM4"),
                                      Mechanotransduction_contractility = c("ABL1", "ACTA1", "ACTC1", "ACTG1", "ACTG2", "ACTN1", "ACTN2", "ACTN3", "ACTN4", "ARF6", "ARHGAP28", "ARHGAP6", "ARHGAP8", "ARHGDIB", "ARHGEF10", "ARHGEF10L", "ARHGEF15", "ARHGEF17", "ARHGEF18", "ARHGEF2", "ARHGEF5", "ARHGEF7", "ARPC4", "ARPC5", "ARPCIA", "CDC42", "CDC42BPA", "CDC42BPB", "CDC42BPG", "CDC42EP2", "CDC42EP5", "CFL1", "FLNA", "FLNB", "KIF14", "LIMCH1", "LIMK1", "MYH10", "MYH11", "MYH13", "MYH14", "MYH15", "MYH2", "MYH3", "MYH4", "MYH7B", "MYH8", "MYH9", "MYL2", "MYL3", "MYL4", "MYL5", "MYL6", "MYL6B", "MYL7", "MYL9", "MYLK2", "MYLK3", "MYLPF", "MYO10", "MYO15A", "MYO16", "MYO18A", "MYO18B", "MYO19", "MYO1A", "MYO1B", "MYO1C", "MYO1D", "MYO1E", "MYO1F", "MYO1G", "MYO1H", "MYO3A", "MYO3B", "MYO5A", "MYO5B", "MYO5C", "MYO6", "MYO7A", "MYO7B", "MYO9A", "MYO9B", "MYOC", "MYOM1", "MYOM2", "MYOZ1", "MYOZ2", "PIEZO1", "RAC1", "RAC2", "RAC3", "RACGAP1", "RHOA", "RHOD", "RHOF", "RHOJ", "ROCK1", "ROCK2", "SRC", "SRF", "TMOD1", "TMOD2", "TMOD3", "TMOD4", "TNNT1", "TNNT2", "TNNT3", "TPM1", "TRIOBP", "TRIP6", "TRPM7", "VIM", "WASF1", "WASF2", "WASF3", "WASL", "WDR1", "WHAMM", "WIPF1")
)

mesen_fa2 <- readRDS("data/mesen_fa2.RDS")
mesen_fa2$time <- ifelse(mesen_fa2$time %in% c('PCW8early','PCW8late'),'PCW8',as.character(mesen_fa2$time)) 
mesen_fa2 <- subset(mesen_fa2,subset = subanno %in% c('Mesen0','Mesen2','Mesen6'))
mesen_fa2 <- subset(mesen_fa2,subset = time %in% c('PCW8','PCW9', 'PCW11', 'PCW13', 'PCW18', 'PCW31'))
DimPlot(mesen_fa2,reduction = 'fa2',group.by = 'subanno')

mesen_fa2 <- AddModuleScore_UCell(mesen_fa2, features = matrix_stiffness_gene_modules, assay="RNA", name = "")

metadf <- mesen_fa2@meta.data
metadf <- metadf[!(metadf$time == 'PCW9'& metadf$subanno == 'Mesen6'),]
metadf <- metadf[!(metadf$time == 'PCW8'& metadf$subanno == 'Mesen6'),]
metadf$id <- row.names(metadf)
metadf$time = factor(metadf$time,levels = c('PCW8','PCW9', 'PCW11', 'PCW13', 'PCW18', 'PCW31'))

modulename = c('ECM_structural_core','ECM_crosslinking_maturation','Mechanotransduction_contractility')

df <- melt(metadf,id.vars = c('id','time','subanno'),measure.vars = modulename,
           value.name = 'ucell')
table(df$subanno,df$time)

result <- aggregate(metadf, 
                    by=list(metadf$time,metadf$subanno), FUN=mean)

result <- result[,c('Group.1','Group.2',modulename)]

dfave <- melt(result,id.vars = c('Group.1','Group.2'),measure.vars = modulename,
              value.name = 'ucell')
colnames(dfave)[1:2] <- c('time','subanno')

ratio_data <- dfave %>%
  pivot_wider(
    id_cols = c(time, variable),
    names_from = subanno,
    values_from = ucell
  ) %>%
  # 过滤掉 Mesen0 或 Mesen2 缺失的时间点
  filter(!is.na(Mesen0) & !is.na(Mesen2)) %>%
  mutate(
    ratio = Mesen2 / Mesen0,
    time = time)


p2 <- ggplot(ratio_data, aes(x =  time, y = ratio,group = 'variable')) + 
  geom_line(size =1,colour = '#9c3587')+ 
  geom_point(colour = '#9c3587',size =3)+
  theme_bw()+
  ylim(0.8,1.5)+
  labs(y = 'Mesen2/Mesen0')+
  facet_wrap(~variable)
p2
ggsave(p2,filename = "figureS6/FigS6b_matrix_stiffness_gene_modules_ratio.pdf",width = 24, height = 3, units = "in", device='pdf')


#figureS6c matrix_stiffness_gene_modules---------------------------------------------------------------------
library(UCell)
matrix_stiffness_gene_modules <- list(ECM_structural_core = c("ACAN", "APCS", "ASPN", "BCAN", "BGN", "C1QL1", "C1QL2", "C1QL3", "C1QL4", "C1QTNF1", "C1QTNF2", "C1QTNF3", "C1QTNF5", "C1QTNF6", "C1QTNF7", "C1QTNF8", "C1QTNF9", "C1QTNF9B", "CCBE1", "CHADL", "CILP", "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", "COL19A1", "COL1A2", "COL20A1", "COL21A1", "COL22A1", "COL23A1", "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1", "COL2A1", "COL3A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A6", "COL5A1", "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6", "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COLEC10", "COLEC11", "COLEC12", "COLGALT1", "COLQ", "COMP", "DCN", "DMP1", "DPT", "ECM1", "ECM2", "EFEMP1", "EFEMP2", "ELN", "EMILIN1", "EMILIN2", "EMILIN3", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FBN3", "FGA", "FGB", "FGG", "FMOD", "FN1", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6", "HAPLN1", "HAPLN2", "HAPLN3", "HAPLN4", "HMCN1", "HMCN2", "HSPG2", "IBSP", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", "LAMB4", "LAMC1", "LAMC2", "LAMC3", "LUM", "MATN1", "MATN2", "MATN3", "MATN4", "MFAP2", "MFAP4", "MFAP5", "MGP", "MMP23B", "NID1", "NID2", "NTN1", "NTN3", "NTN4", "OGN", "POSTN", "PRELP", "SFTPA", "SFTPA2", "SFTPD", "SLIT2", "SMOC1", "SMOC2", "SPARC", "SPON1", "SPP1", "SPP2", "SRPX", "SRPX2", "THBS1", "THBS2", "THBS3", "THBS4", "TINAG", "TINAGL1", "TNC", "TNN", "TNR", "TNXA", "TNXB", "VIT", "VTN", "VWA1", "VWA2", "VWF"),
                                      ECM_crosslinking_maturation = c("CRTAP", "FKBP10", "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "P3H2", "P3H3", "P4HA1", "P4HA2", "P4HA3", "P4HB", "PCOLCE", "PCOLCE2", "PLOD3", "PXDN", "SERPINH1", "TGM2", "TGM4"),
                                      Mechanotransduction_contractility = c("ABL1", "ACTA1", "ACTC1", "ACTG1", "ACTG2", "ACTN1", "ACTN2", "ACTN3", "ACTN4", "ARF6", "ARHGAP28", "ARHGAP6", "ARHGAP8", "ARHGDIB", "ARHGEF10", "ARHGEF10L", "ARHGEF15", "ARHGEF17", "ARHGEF18", "ARHGEF2", "ARHGEF5", "ARHGEF7", "ARPC4", "ARPC5", "ARPCIA", "CDC42", "CDC42BPA", "CDC42BPB", "CDC42BPG", "CDC42EP2", "CDC42EP5", "CFL1", "FLNA", "FLNB", "KIF14", "LIMCH1", "LIMK1", "MYH10", "MYH11", "MYH13", "MYH14", "MYH15", "MYH2", "MYH3", "MYH4", "MYH7B", "MYH8", "MYH9", "MYL2", "MYL3", "MYL4", "MYL5", "MYL6", "MYL6B", "MYL7", "MYL9", "MYLK2", "MYLK3", "MYLPF", "MYO10", "MYO15A", "MYO16", "MYO18A", "MYO18B", "MYO19", "MYO1A", "MYO1B", "MYO1C", "MYO1D", "MYO1E", "MYO1F", "MYO1G", "MYO1H", "MYO3A", "MYO3B", "MYO5A", "MYO5B", "MYO5C", "MYO6", "MYO7A", "MYO7B", "MYO9A", "MYO9B", "MYOC", "MYOM1", "MYOM2", "MYOZ1", "MYOZ2", "PIEZO1", "RAC1", "RAC2", "RAC3", "RACGAP1", "RHOA", "RHOD", "RHOF", "RHOJ", "ROCK1", "ROCK2", "SRC", "SRF", "TMOD1", "TMOD2", "TMOD3", "TMOD4", "TNNT1", "TNNT2", "TNNT3", "TPM1", "TRIOBP", "TRIP6", "TRPM7", "VIM", "WASF1", "WASF2", "WASF3", "WASL", "WDR1", "WHAMM", "WIPF1")
)

load('data/ut_15_008.Rdata')
object <- AddModuleScore_UCell(object, features = matrix_stiffness_gene_modules,name = '')

coords <- GetTissueCoordinates(object)[,1:2]
colnames(coords)<-c("FDG1","FDG2")
metadf <- object@meta.data
myeloid <- merge(coords,metadf,by = 0)

myeloid <- myeloid[myeloid$FDG1 > -30000 & myeloid$FDG1 < -15000 & myeloid$FDG2 < 24000 & myeloid$FDG2 > 16000,]

pdf('figureS6/figS6c_MP_ut_15_matrix_stiffness_gene_modules.pdf', width = 7, height = 5)
p1<-ggplot(myeloid)+
  geom_point(aes_string(x = 'FDG1', y = 'FDG2',fill='ECM_structural_core'),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  theme_bw()+
  scale_fill_distiller(palette = "RdYlGn", direction = -1,limits = c(0, 0.23),oob = scales::oob_squish)+
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

p1<-ggplot(myeloid)+
  geom_point(aes_string(x = 'FDG1', y = 'FDG2',fill='ECM_crosslinking_maturation'),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  theme_bw()+
  scale_fill_distiller(palette = "RdYlGn", direction = -1,limits = c(0, 0.12),oob = scales::oob_squish)+
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


p1<-ggplot(myeloid)+
  geom_point(aes_string(x = 'FDG1', y = 'FDG2',fill='Mechanotransduction_contractility'),shape = 21,stroke = 0.2,size = .6,color = rgb(0,0,0,0))+
  theme_bw()+
  scale_fill_distiller(palette = "RdYlGn", direction = -1,limits = c(0, 0.15),oob = scales::oob_squish)+
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


dev.off()

#figureS6d MP matrix_stiffness_gene_modules----------------------------------------------------------------------------
library(pheatmap)
library(tidyr)

enrichlist <- read.csv('data/MP_matrix_stiffness_gene_modules_NMF50.csv')

enrichlist$MP <- paste0('MP',enrichlist$MP)
enrichlist <- enrichlist[,c('Description','pvalue','MP')]
enrich_width <- spread(enrichlist, key = "MP",
                       value = "pvalue")

row.names(enrich_width) <- enrich_width$Description
enrich_width <- enrich_width[,-1]
enrich_width <- enrich_width[unique(enrichlist$Description),unique(enrichlist$MP)]

enrich_width[is.na(enrich_width)] <- 1
enrich_width <- -log10(enrich_width)


p1 <- pheatmap(enrich_width,cluster_rows = F,cluster_cols = F,cellwidth = 10,
               show_rownames = T,show_colnames = T,
               color = colorRampPalette(colors = c( "white", "red"))(100),
               treeheight_row = 0,treeheight_col = 0,
               breaks = seq(0, 2, length.out = 101)) 
p1


pdf('figureS6/figS6d_matrix_stiffness_gene_modules_enrichment.pdf', width = 7, height = 2)
p1
dev.off()

#figureS6e matrix_stiffness_gene_modules---------------------------------------------------------------------
library(UCell)
matrix_stiffness_gene_modules <- list(ECM_structural_core = c("ACAN", "APCS", "ASPN", "BCAN", "BGN", "C1QL1", "C1QL2", "C1QL3", "C1QL4", "C1QTNF1", "C1QTNF2", "C1QTNF3", "C1QTNF5", "C1QTNF6", "C1QTNF7", "C1QTNF8", "C1QTNF9", "C1QTNF9B", "CCBE1", "CHADL", "CILP", "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", "COL19A1", "COL1A2", "COL20A1", "COL21A1", "COL22A1", "COL23A1", "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1", "COL2A1", "COL3A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A6", "COL5A1", "COL5A2", "COL5A3", "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6", "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COLEC10", "COLEC11", "COLEC12", "COLGALT1", "COLQ", "COMP", "DCN", "DMP1", "DPT", "ECM1", "ECM2", "EFEMP1", "EFEMP2", "ELN", "EMILIN1", "EMILIN2", "EMILIN3", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FBN3", "FGA", "FGB", "FGG", "FMOD", "FN1", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6", "HAPLN1", "HAPLN2", "HAPLN3", "HAPLN4", "HMCN1", "HMCN2", "HSPG2", "IBSP", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", "LAMB4", "LAMC1", "LAMC2", "LAMC3", "LUM", "MATN1", "MATN2", "MATN3", "MATN4", "MFAP2", "MFAP4", "MFAP5", "MGP", "MMP23B", "NID1", "NID2", "NTN1", "NTN3", "NTN4", "OGN", "POSTN", "PRELP", "SFTPA", "SFTPA2", "SFTPD", "SLIT2", "SMOC1", "SMOC2", "SPARC", "SPON1", "SPP1", "SPP2", "SRPX", "SRPX2", "THBS1", "THBS2", "THBS3", "THBS4", "TINAG", "TINAGL1", "TNC", "TNN", "TNR", "TNXA", "TNXB", "VIT", "VTN", "VWA1", "VWA2", "VWF"),
                                      ECM_crosslinking_maturation = c("CRTAP", "FKBP10", "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "P3H2", "P3H3", "P4HA1", "P4HA2", "P4HA3", "P4HB", "PCOLCE", "PCOLCE2", "PLOD3", "PXDN", "SERPINH1", "TGM2", "TGM4"),
                                      Mechanotransduction_contractility = c("ABL1", "ACTA1", "ACTC1", "ACTG1", "ACTG2", "ACTN1", "ACTN2", "ACTN3", "ACTN4", "ARF6", "ARHGAP28", "ARHGAP6", "ARHGAP8", "ARHGDIB", "ARHGEF10", "ARHGEF10L", "ARHGEF15", "ARHGEF17", "ARHGEF18", "ARHGEF2", "ARHGEF5", "ARHGEF7", "ARPC4", "ARPC5", "ARPCIA", "CDC42", "CDC42BPA", "CDC42BPB", "CDC42BPG", "CDC42EP2", "CDC42EP5", "CFL1", "FLNA", "FLNB", "KIF14", "LIMCH1", "LIMK1", "MYH10", "MYH11", "MYH13", "MYH14", "MYH15", "MYH2", "MYH3", "MYH4", "MYH7B", "MYH8", "MYH9", "MYL2", "MYL3", "MYL4", "MYL5", "MYL6", "MYL6B", "MYL7", "MYL9", "MYLK2", "MYLK3", "MYLPF", "MYO10", "MYO15A", "MYO16", "MYO18A", "MYO18B", "MYO19", "MYO1A", "MYO1B", "MYO1C", "MYO1D", "MYO1E", "MYO1F", "MYO1G", "MYO1H", "MYO3A", "MYO3B", "MYO5A", "MYO5B", "MYO5C", "MYO6", "MYO7A", "MYO7B", "MYO9A", "MYO9B", "MYOC", "MYOM1", "MYOM2", "MYOZ1", "MYOZ2", "PIEZO1", "RAC1", "RAC2", "RAC3", "RACGAP1", "RHOA", "RHOD", "RHOF", "RHOJ", "ROCK1", "ROCK2", "SRC", "SRF", "TMOD1", "TMOD2", "TMOD3", "TMOD4", "TNNT1", "TNNT2", "TNNT3", "TPM1", "TRIOBP", "TRIP6", "TRPM7", "VIM", "WASF1", "WASF2", "WASF3", "WASL", "WDR1", "WHAMM", "WIPF1")
)


sp_fa2 <- readRDS('data/sp_fa2.RDS')
sp_fa2 <- AddModuleScore_UCell(sp_fa2, features = matrix_stiffness_gene_modules,name = '')

myeloid_fdg_coor <- sp_fa2@reductions[["fa2"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- sp_fa2@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)

pdf('figureS6/figS6e_MP_fa2_matrix_stiffness_gene_modules.pdf', width = 8, height = 6)
for (i in c('ECM_structural_core','ECM_crosslinking_maturation','Mechanotransduction_contractility')) {
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
    )
  
  print(p1)
}
dev.off()



#figureS6h snATAC tsne----------------------------------------------------------------
library(ggplot2)
library(Seurat)

load('data/ATAC_anno.Rdata')
atac <- RunTSNE(object = pbmc, reduction = 'lsi', dims = 2:30)

pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff",'#fefee3','#d68c45', '#ffc9b9','#ac6f82'
)

names(pal11) <-  c("Endo","Epi.MD", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen",'Mesen0','Mesen2','Mesen6','Others')
myeloid_fdg_coor <- atac@reductions[["tsne"]]@cell.embeddings
colnames(myeloid_fdg_coor)<-c("FDG1","FDG2")
metadf <- atac@meta.data
myeloid <- merge(myeloid_fdg_coor,metadf,by = 0)
myeloid$plot_order <- ifelse(myeloid$predicted.id == "Mesen2", 2, 
                             ifelse(myeloid$predicted.id == "Mesen0", 1, 0))
myeloid <- myeloid[order(myeloid$plot_order),]

p1<-ggplot(myeloid)+
  geom_point(aes(x = FDG1, y = FDG2,fill=predicted.id),shape = 21,stroke = 0.2,size = 3)+
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
ggsave(p1,filename = "figureS6/figS6h_atac_Nolegend.pdf",width = 6, height = 6, units = "in", device='pdf')
