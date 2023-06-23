#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggrepel)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_exp_vs_KD <- function(df, graphname){
  textsize <- 7
  print (paste('correlation for:', graphname, cor(df$log10_Kd, df$Exp_score)))
  class_levels <- c('CR9114 mutant', 'Nonsense', 'Influenza (HA stem)', 'Influenza (non-HA stem)', 'Non-influenza', 'WT')
  df <- df %>%
          mutate(class=factor(class, levels=class_levels))
  palette <- c(qualpal(n = 8, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  palette <- c(palette[1:3], palette[7], 'gray', palette[8])
  p <-  ggplot(df, aes(x=Exp_score, y=log10_Kd, color=class, size=log10(avg_input_freq))) +
          geom_point(alpha=0.7, pch=16) +
          scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-5.5, -2.5),
                                range = c(0.1, 0.8),
                                breaks = c(-5, -4, -3),
                                labels = c(expression(bold('10'^'-5')), expression(bold('10'^'-4')), expression(bold('10'^'-3')))) +
          scale_color_manual('',values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab(expression(bold(log['10']~K['D']))) +
          xlab("Expression score")
  ggsave(graphname,p,width=3.2, height=1.9, bg='white', dpi=600)
  }

assign_epi <- function(class, epi){
  if (class=='nonsense'){return ('Nonsense')}
  else if (class=='CR9114 mutant'){return ('CR9114 mutant')}
  else if (class=='WT'){return ('WT')}
  else if (epi=='HA stem'){return ('Influenza (HA stem)')}
  else if (class=='influenza'){return ('Influenza (non-HA stem)')}
  else {return ('Non-influenza')}
  }


df_exp <- read_csv('result/CDRH3_express_table_summary.csv')
df_bind <- read_csv('result/CDRH3_KD_table_summary_class.csv')
df <- inner_join(df_bind,df_exp) %>%
  mutate(log10_Kd=log10(Kd)) %>%
  filter((log10_Kd < -8 & p.value < 0.2) | (log10_Kd >= -8)) %>%
  mutate(class = mapply(assign_epi, class, epi))
plot_exp_vs_KD(df, 'graph/CDRH3_KD_vs_exp.png')
