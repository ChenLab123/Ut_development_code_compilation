library(readxl)
library(tidyr)
library(viridis)


#DES MP--------------------------------------------------------------------------------
filelist <- list.files('data/实验数据/QUPATH OUTPUT/DES MP')

DESline <- data.frame()

for (i in filelist) {
  df <- read.table(paste0('data/实验数据/QUPATH OUTPUT/DES MP/',i,'/ANNOTATION.txt'),header = T,sep = '\t')
  df$time <- i
  max_len <- max(nchar(names(df)))
  longest_columns <- names(df)[nchar(names(df)) == max_len]
  df <- df[,c('time','Name','Num.Detections',longest_columns)]
  colnames(df) <- c('time','Name','Num.Detections','FourGenes')
  DESline <- rbind(DESline,df)
}

DESline <- DESline[DESline$Name %in% c('LEFT','RIGHT'),]
DESline$percent <- DESline$FourGenes/DESline$Num.Detections
p1 <- ggplot(DESline, aes(x =  time, y = percent,group = Name,colour = Name,fill = Name)) + 
  geom_line(size =1)+ 
  geom_point(size =3)+
  facet_wrap(~Name)+
  theme_bw()+
  labs(x = 'FourGenes percent')


filelist <- list.files('data/实验数据/QUPATH OUTPUT/DES MP')

DESbar <- data.frame()

for (i in filelist) {
  df <- read.table(paste0('data/实验数据/QUPATH OUTPUT/DES MP/',i,'/DETECTION.txt'),header = T,sep = '\t')
  df$time <- i
  colnames(df) <- gsub('DES.','DESMIN.',colnames(df), fixed = TRUE)
  DESbar <- rbind(DESbar,df)
}
DESbar <- DESbar[DESbar$Parent %in% c('LEFT','RIGHT'),]
colnames(DESbar)
DESbar <- DESbar[,c('time','Parent','Cell..520.ACTG2.mean','Cell..570.CNN1.mean',
                    'Cell..690.NCAM1.mean','Cell..620.DESMIN.mean')]

DESbar$time <- factor(DESbar$time,levels = unique(DESbar$time))

DESbar <- gather(DESbar, key = "Gene", value = "value",
                 -'time', -'Parent')

df_stats <- DESbar %>%
  group_by(Gene, time, Parent) %>%
  summarise(
    sd = sd(value),
    mean = mean(value),
    n = n(),
    se = sd / sqrt(n),  # 这里用的 sd 和 n 是这一行之前定义的
    .groups = "drop"
  )


p2 <- ggplot(df_stats, aes(x = time, y = mean, fill = Gene)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.9),
                width = 0.2, size = 0.8) +
  facet_wrap(~Parent) +
  labs(y = "Mean ± SE")+
  theme_classic()+
  scale_fill_viridis(discrete = T)
p2

pdf('实验结果/DES_plot.pdf', width = 12, height = 8)
p1/p2+ plot_layout(heights = c(1,5))
dev.off()

#TRANS--------------------------------------------------------------------------------
filelist <- list.files('data/实验数据/QUPATH OUTPUT/TRANS')

TRANSline <- data.frame()

for (i in filelist) {
  df <- read.table(paste0('data/实验数据/QUPATH OUTPUT/TRANS/',i,'/ANNOTATION.txt'),header = T,sep = '\t')
  df$time <- i
  max_len <- max(nchar(names(df)))
  longest_columns <- names(df)[nchar(names(df)) == max_len]
  print(longest_columns)
  df <- df[,c('time','Name','Num.Detections',longest_columns)]
  colnames(df) <- c('time','Name','Num.Detections','FourGenes')
  TRANSline <- rbind(TRANSline,df)
}

TRANSline <- TRANSline[TRANSline$Name %in% c('LEFT','RIGHT'),]
TRANSline$percent <- TRANSline$FourGenes/TRANSline$Num.Detections
p1 <- ggplot(TRANSline, aes(x =  time, y = percent,group = Name,colour = Name,fill = Name)) + 
  geom_line(size =1)+ 
  geom_point(size =3)+
  facet_wrap(~Name)+
  theme_bw()+
  labs(x = 'FourGenes percent')
p1

filelist <- list.files('data/实验数据/QUPATH OUTPUT/TRANS')

TRANSbar <- data.frame()

for (i in filelist) {
  df <- read.table(paste0('data/实验数据/QUPATH OUTPUT/TRANS/',i,'/DETECTION.txt'),header = T,sep = '\t')
  df$time <- i
  TRANSbar <- rbind(TRANSbar,df)
}
TRANSbar <- TRANSbar[TRANSbar$Parent %in% c('LEFT','RIGHT'),]
colnames(TRANSbar)
TRANSbar <- TRANSbar[,c('time','Parent','Cell..520.LUMICAN.mean','Cell..570.CXCL12.mean',
                        'Cell..690.SFRP1.mean','Cell..620.IGFBP3.mean')]

TRANSbar$time <- factor(TRANSbar$time,levels = unique(TRANSbar$time))

TRANSbar <- gather(TRANSbar, key = "Gene", value = "value",
                   -'time', -'Parent')

df_stats <- TRANSbar %>%
  group_by(Gene, time, Parent) %>%
  summarise(
    sd = sd(value),
    mean = mean(value),
    n = n(),
    se = sd / sqrt(n),  # 这里用的 sd 和 n 是这一行之前定义的
    .groups = "drop"
  )


p2 <- ggplot(df_stats, aes(x = time, y = mean, fill = Gene)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.9),
                width = 0.2, size = 0.8) +
  facet_wrap(~Parent) +
  labs(y = "Mean ± SE")+
  theme_classic()+
  scale_fill_brewer(palette = "Set1") 
p2


pdf('实验结果/TRANS_plot.pdf', width = 12, height = 8)
p1/p2+ plot_layout(heights = c(1,5))
dev.off()

