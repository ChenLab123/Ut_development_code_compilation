library(ggplot2)
library(Seurat)

load('data/mus_scRNA.Rdata')
pal11 <- c(
  "#debcff", "#caffbf", "#96e8ff","#ffadad", 
  "#a0c4ff", "#F0E442","#ffd6a5",  "#ffc6ff", "#ffa7dc","#74B346",
  "#00c2ff",'#4c956c'
)

names(pal11) <-  c("Endo","Epi", "Eryth" ,"ESC", "ESR1.Mesen", "Immu","Mesen","NCC",
                   "prolif","SMC","THY1.Mesen",'VSMC')

#figureS7 tsne-----------------------------------------------------------------------
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
  )

p1
ggsave(p1,filename = "figureS7/figS7a_scRNA_withlegend.pdf",width = 6, height = 6, units = "in", device='pdf')

#figS7b alluvial_plot-----------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggalluvial)

load('data/mus_scRNA.Rdata')
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


pdf('figureS7/figS7b_alluvial_plot.pdf', width = 10, height = 6)
p1
dev.off()


#fig5 circle-----------------------------------------------------------------------
scRNA@meta.data[["time1"]] <- ifelse(scRNA@meta.data[["time"]] %in% c('PND1','PND5','PND10','PND15'),'PND','Adu')

df <- as.data.frame(table(scRNA$mus_anno,scRNA$time1)) 

colnames(df) <- c('type','category','count')


data <- df[df$type == 'ESC',]

data$fraction <- data$count / sum(data$count)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

data$labelPosition <- (data$ymax + data$ymin) / 2

data$label <- paste0(round(data$fraction,2))

p1 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  scale_fill_brewer(palette=2) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  ggtitle('ESC')

data <- df[df$type == 'Mesen',]

data$fraction <- data$count / sum(data$count)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

data$labelPosition <- (data$ymax + data$ymin) / 2

data$label <- paste0(round(data$fraction,2))

p2 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  scale_fill_brewer(palette=2) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  ggtitle('Mesen')

data <- df[df$type == 'SMC',]

data$fraction <- data$count / sum(data$count)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

data$labelPosition <- (data$ymax + data$ymin) / 2

data$label <- paste0(round(data$fraction,2))

p3 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  scale_fill_brewer(palette=2) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  ggtitle('SMC')

data <- df[df$type == 'VSMC',]

data$fraction <- data$count / sum(data$count)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

data$labelPosition <- (data$ymax + data$ymin) / 2

data$label <- paste0(round(data$fraction,2))

p4 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  scale_fill_brewer(palette=2) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "right")+
  ggtitle('VSMC')

(p1+p2)/(p3+p4)

pdf('figureS7/figS7c_donut.pdf', width = 5, height = 5)
(p1+p2)/(p3+p4)
dev.off()

