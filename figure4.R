library(ggplot2)
#figure4b matrix_stiffness_gene_modules----------------------------------------------------------------------------
library(Seurat)
library(UCell)
library(dplyr)
library(reshape2)
library(tidyr)
library(patchwork)

col.mesen <- c('#e6d8a0','#d68c45', '#ffc9b9')
names(col.mesen) <- c('Mesen0','Mesen2','Mesen6')

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


p1 <- ggplot(dfave, aes(x =  time, y = ucell,group = subanno,colour = subanno,fill = subanno)) + 
  geom_line(size =1)+ 
  geom_point(size =3)+
  scale_fill_manual(values = as.character(col.mesen),breaks = names(col.mesen))+
  scale_color_manual(values = as.character(col.mesen),breaks = names(col.mesen))+
  theme_bw()+
  facet_wrap(~variable, scales = "free_y")
p1

ggsave(p1,filename = "figure4/Fig4b_ECM_matrix_stiffness_gene_modules.pdf",width = 24, height = 4, units = "in", device='pdf')

#figure4e dorothea----------------------------------------------------------------------------
library(tibble)
library(pheatmap)
load('data/dorothea.Rdata')
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
  top_n(114, var) %>%
  distinct(tf)
highly_variable_tfs <- rbind(highly_variable_tfs,'TWIST1')

summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

summarized_viper_scores_df[1:4,1:4]
summarized_viper_scores_df <- t(summarized_viper_scores_df)
summarized_viper_scores_df <- summarized_viper_scores_df[,c('VSMC','Mesen2','Mesen0','Mesen6','DES.SMC','ESC')]

palette_length = 100
my_color = colorRampPalette(c("#3e5aa8", "white","#ff6361"))(100)

p1 <- pheatmap(summarized_viper_scores_df,fontsize=14, 
                       fontsize_row = 10, cluster_cols = F,
                       color=my_color, breaks = seq(-0.5, 0.5, length.out = 101), 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

pdf('figure4/fig4e_dorothea.pdf', width = 8, height = 6)
p1
dev.off()

#figure4h lollipop----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(hrbrthemes)

data <- read.csv("data/转录效率_mesen6.csv")

# Reorder the data
data <- data %>%
  arrange(avg_expr) %>%
  mutate(x=factor(tf,tf))

# Plot
p1 <- ggplot(data, aes(x=tf, y=avg_expr)) +
  geom_segment(
    aes(x=x, xend=x, y=0, yend=avg_expr), 
    color=ifelse(data$x %in% c("TWIST1"), "orange", "grey"), 
    size=ifelse(data$x %in% c("TWIST1"), 1.3, 1.3)
  ) +
  geom_point(
    color=ifelse(data$x %in% c("TWIST1"), "orange", "grey"), 
    size=ifelse(data$x %in% c("TWIST1"), 5, 2)
  ) +
  theme_classic() +
  coord_flip() +
  theme(
    legend.position="none"
  ) +
  xlab("") +
  ylab("avg_expr") +
  ggtitle("Mesen6")
print(p1)


pdf('figure4/fig4h_lollipop.pdf', width = 5, height = 9)
p1
dev.off()

#figure4i ipsc----------------------------------------------------------------------------
esdf <- read.csv('data/ATAC_esdf.csv')

color_map <- setNames(
  rep("grey88", length(unique(esdf$Type))),
  unique(esdf$Type)
)
color_map["Smooth.muscle.cells"] <- "red"


p1 <- ggplot(esdf, aes(x = num, y = enrichment_scores, color = Type)) + 
  geom_line(alpha = 0.8) +
  scale_color_manual(values = color_map) +
  theme_classic() +
  theme(legend.position = "none")+
  labs(title = 'Smooth.muscle.cells')
print(p1)


pdf('figure4/fig4i_ipsc_SMC.pdf', width = 5, height = 4)
p1
dev.off()


#figure4f ATAC peaks----------------------------------------------------------------------------
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(scCustomize)

pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff",'#fefee3','#d68c45', '#ffc9b9','#ac6f82'
)

names(pal11) <-  c("Endo","Epi.MD", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen",'Mesen0','Mesen2','Mesen6','Others')

load("data/ATAC_chrom.Rdata")
DefaultAssay(atac) <- 'peaks'

p1 <- CoveragePlot(
  object = atac,
  region = 'ACTA2',
  extend.upstream = 20000,
  extend.downstream = 20000, 
  ymax =15,
  group.by ='predicted.id'
)  
p1 <- p1 & scale_fill_manual(values = as.character(pal11),breaks = names(pal11))
p1

ggsave(p1,filename = "figure4/fig4f_ACTA2_peaks.pdf",width = 6, height = 4, units = "in", device='pdf')

#figure4g motifs----------------------------------------------------------------------------
p1 <- MotifPlot(
  object = atac,
  motifs = 'MA1123.2',
  assay = 'peaks'
)
ggsave(p1,filename = "figure4/fig4g_TWIST1_motif.pdf",width = 4, height = 2, units = "in", device='pdf')

p1 <-MotifPlot(
  object = atac,
  motifs = 'MA0083.3',
  assay = 'peaks'
)
ggsave(p1,filename = "figure4/fig4g_SRF_motif.pdf",width = 4, height = 2, units = "in", device='pdf')

#figure4k TWIST1 peaks----------------------------------------------------------------------------
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(scCustomize)

pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff",'#fefee3','#d68c45', '#ffc9b9','#ac6f82'
)

names(pal11) <-  c("Endo","Epi.MD", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen",'Mesen0','Mesen2','Mesen6','Others')

load("data/ATAC_chrom.Rdata")
DefaultAssay(atac) <- 'peaks'

p1 <- CoveragePlot(
  object = atac,
  region = 'TWIST1',
  extend.upstream = 20000,
  extend.downstream = 20000, 
  ymax =20,
  group.by ='predicted.id'
)  
p1 <- p1 & scale_fill_manual(values = as.character(pal11),breaks = names(pal11))
p1
ggsave(p1,filename = "figure4/fig4k_TWIST1_peaks.pdf",width = 6, height = 4, units = "in", device='pdf')

#figure4l motifs----------------------------------------------------------------------------
p1 <- MotifPlot(
  object = atac,
  motifs = 'MA0090.3',
  assay = 'peaks'
)
ggsave(p1,filename = "figure4/fig4l_TEAD1_motif.pdf",width = 4, height = 2, units = "in", device='pdf')

p1 <- MotifPlot(
  object = atac,
  motifs = 'MA1121.1',
  assay = 'peaks'
)
ggsave(p1,filename = "figure4/fig4l_TEAD2_motif.pdf",width = 4, height = 2, units = "in", device='pdf')

p1 <- MotifPlot(
  object = atac,
  motifs = 'MA0808.1',
  assay = 'peaks'
)
ggsave(p1,filename = "figure4/fig4l_TEAD3_motif.pdf",width = 4, height = 2, units = "in", device='pdf')

p1 <-MotifPlot(
  object = atac,
  motifs = 'MA0809.2',
  assay = 'peaks'
)
ggsave(p1,filename = "figure4/fig4l_TEAD4_motif.pdf",width = 4, height = 2, units = "in", device='pdf')
