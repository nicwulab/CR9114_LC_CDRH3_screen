#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(stringr)
require(cowplot)

plot_fitness_heatmap <- function(df, WTresibox, start_resi, end_resi){
  textsize <- 6
  df <- df %>%
                     filter(Pos >= start_resi & Pos <= end_resi)
  WTresibox     <- WTresibox %>%
                     filter(Pos >= start_resi & Pos <= end_resi) %>%
                     mutate(x=x-min(x)+1)
  p <-  ggplot() +
          geom_tile(data=df,aes(x=resi,y=aa,fill=log10_Kd)) +
          scale_fill_gradientn('KD',
                colours=c("red","white","white","blue","blue"),
                limits=c(-10,-4.9),
                values=rescale(c(-10,-9.3,-9.1,-8,-4.9)),
                breaks=c(-10,-9,-8,-7,-6,-5),
                labels=c(expression(bold('10'^'-10')), expression(bold('10'^'-9')), expression(bold('10'^'-8')), 
                         expression(bold('10'^'-7')), expression(bold('10'^'-6')), expression(bold('10'^'-5'))),
                guide="colorbar",
                na.value="grey") +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=7,face="bold"),
                axis.line = element_line(colour = 'black', linewidth = 0),
                legend.text = element_text(hjust = -1),
                panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
          guides(fill = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                       label.theme=element_text(size=7,face="bold",colour='black'),
                                       frame.colour="black",
                                       frame.linewidth = 1,
                                       ticks = TRUE,
                                       ticks.colour = "black",
                                       barwidth = 0.5, barheight = 6, title=expression(bold(log['10']~K['d'])))) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.2) +
          scale_x_discrete(labels=c('A93','R94','H95','G96','N97','Y98','Y99','Y100','Y100a','S100b','G100c','M100d','D101','V102')) +
          xlab("") +
          ylab("amino acid")
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))

df <- read_csv('result/CDRH3_KD_table_summary.csv') %>%
  filter(grepl('CR9114',ID)) %>%
  mutate(log10_Kd=log10(Kd)) %>%
  filter((log10_Kd < -8 & p.value < 0.2) | (log10_Kd >= -8)) %>%
  mutate(Mutation=gsub('CR9114_',"",ID)) %>%
  filter(Mutation != 'WT') %>%
  mutate(resi=str_sub(Mutation,1,-2)) %>%
  mutate(aa=str_sub(Mutation,-1,-1)) %>%
  filter(aa %in% aa_level) %>%
  mutate(aa=factor(aa,levels=aa_level)) %>%
  complete(resi, aa) %>%
  mutate(Pos=str_sub(resi,2,-1)) %>%
  mutate(Pos=as.numeric(as.character(Pos))) %>%
  arrange(Pos) %>%
  mutate(resi=factor(resi,levels=unique(resi))) %>%
  mutate(log10_Kd=case_when(str_sub(resi,1,1)==aa ~ log10(5.19e-10), TRUE ~ log10_Kd)) %>%
  mutate(Mutation=paste(resi,aa,sep='')) %>%
  select(Mutation, resi, Pos, aa, log10_Kd)

WTresibox  <- df %>%
  select(resi,Pos) %>%
  unique() %>%
  mutate(WT_resi=str_sub(resi,1,1)) %>%
  mutate(x=seq(1,14)) %>%
  mutate(y=match(WT_resi,aa_level)) %>%
  select(resi,WT_resi,Pos,x, y)

print (range(df$log10_Kd,na.rm=TRUE))

p1 <- plot_fitness_heatmap(df, WTresibox, 1, 14)
p <- grid.arrange(p1, nrow=1)
ggsave('graph/CR9114_CDRH3_KD_heatmap.png',p,width=2.2, height=2.5, dpi=600)
