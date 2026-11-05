library(tidyverse)
library(patchwork)
library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)


load("data/JoinLayers_scRNA.rdata")

#figureS1b-----------------------------------------------------------------
qcdf <- combined@meta.data

qcdf <- qcdf[,c("nCount_RNA","nFeature_RNA","time")]

sample_size = qcdf %>% group_by(time) %>% summarize(num=n())
sample_size$myaxis = paste0(sample_size$time, "\n", "n=", sample_size$num)

p1 <- ggplot(qcdf, aes(x=time, y=nCount_RNA, fill=time)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2,outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("") +
  xlab("") +
  coord_cartesian(ylim = c(0, 20000))+ 
  scale_x_discrete(labels = sample_size$myaxis)

p2 <- ggplot(qcdf, aes(x=time, y=nFeature_RNA, fill=time)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2,outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("") +
  xlab("") +
  coord_cartesian(ylim = c(0, 5000))+ 
  scale_x_discrete(labels = sample_size$myaxis)

p1/p2

ggsave(p1/p2,filename = "figureS1/figS1b_scRNA_QC.pdf",width = 10, height = 6, units = "in", device='pdf')
