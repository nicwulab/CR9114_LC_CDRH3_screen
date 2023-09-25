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

plot_dist <- function(df, parameter, graphname, ylab){
  textsize <- 7
  class_levels <- c('CR9114\nmutant', 'Nonsense', 'Influenza\n(HA stem)', 'Influenza\n(non-HA stem)', 'Non-influenza')
  print (graphname)
  for (class_level in class_levels){
    t <- t.test(filter(df,class=='CR9114\nmutant')$param,filter(df,class==class_level)$param)
    print (paste('CR9114\nmutant', class_level, t$p.value))
    }
  palette <- c(qualpal(n = 8, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  palette <- c(palette[1:3], palette[7], 'gray', palette[8])
  WT_log10_Kd <- filter(df, class=='WT')$param
  df <- filter(df, class != 'WT') %>%
          mutate(class=factor(class, levels=class_levels))
  p <-  ggplot() +
          geom_violin(data=df,aes(x=class,y=param),width=1, size=0.3) +
          geom_sina(data=df,aes(x=class,y=param,color=class),
                    pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, alpha=0.3) +
          geom_boxplot(data=df,aes(x=class,y=param),width=0.1, size=0.3, alpha=0, outlier.shape=NA) +
          scale_color_manual('',values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=7,face="bold"),
                legend.position = "none") +
          xlab("") +
          ylab(ylab)
  if (parameter=='KD'){
     print (WT_log10_Kd)
     p <- p+geom_hline(yintercept=WT_log10_Kd, linetype = "dashed")
     }
  ggsave(graphname,p,width=2.8, height=1.8, dpi=600, bg='white')
  }

plot_motif_analysis <- function(df, graphname){
  textsize <- 7
  class_levels <- c('Influenza\n(HA stem)', 'Influenza\n(non-HA stem)', 'Non-influenza')
  df <- filter(df, class %in% class_levels) %>%
          mutate(class=factor(class, levels=class_levels))
  for (class_level in class_levels){
    t <- t.test(filter(df, class==class_level & param=='with Tyr')$log10_Kd, filter(df, class==class_level & param=='without Tyr')$log10_Kd)
    print (paste(graphname, class_level, t$p.value))
    }
  palette <- brewer.pal(3,"Set1")
  WT_log10_Kd <- filter(df, class=='WT')$param
  dodge_value <- 0.8
  p <-  ggplot(df, aes(x=class, y=log10_Kd, color=param, fill=param)) +
          geom_sina(alpha=0.8, size=0.1, pch=16, position=position_dodge(dodge_value), maxwidth=0.4) +
          geom_boxplot(width=0.5, alpha=0.4, size=0.3, outlier.shape=NA, color='black', position=position_dodge(dodge_value)) +
          scale_color_manual('',values=palette,drop=FALSE) +
          scale_fill_manual('',values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=7,face="bold"),
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.height=unit(0.8,"line"),
                legend.key.size = unit(0.4,"line"),
                legend.position = "right") +
          xlab("") +
          ylab('log10 KD')
  ggsave(graphname,p,width=2.8, height=2, dpi=600, bg='white')
  }

assign_epi <- function(class, epi){
  if (class=='nonsense'){return ('Nonsense')}
  else if (class=='CR9114 mutant'){return ('CR9114\nmutant')}
  else if (class=='WT'){return ('WT')}
  else if (epi=='HA stem'){return ('Influenza\n(HA stem)')}
  else if (class=='influenza'){return ('Influenza\n(non-HA stem)')}
  else {return ('Non-influenza')}
  }

analyze_motif <- function(motif){
  if (grepl('Y',motif)){return ('with Tyr')}
  else {return ('without Tyr')}
  }

df_exp <- read_csv('result/CDRH3_express_table_summary.csv')
df_bind <- read_csv('result/CDRH3_KD_table_summary_class.csv')
df <- inner_join(df_bind,df_exp) %>%
  mutate(log10_Kd=log10(Kd)) %>%
  filter((log10_Kd < -8 & p.value < 0.2) | (log10_Kd >= -8)) %>%
  mutate(class = mapply(assign_epi, class, epi)) %>%
  mutate(motif_98 = mapply(analyze_motif, resi98)) %>%
  mutate(motif_CDRH3 = mapply(analyze_motif, seq))

print (range(df$log10_Kd,na.rm=TRUE))
print (length(df$log10_Kd))
print (table(df$class))
print (table(df$epi))
print (mean(filter(df, class=="CR9114\nmutant")$log10_Kd))

df_plot <- rename(df, param=log10_Kd)
plot_dist(df_plot, 'KD', 'graph/CDRH3_distrib_KD.png', expression(bold(log['10']~K['D'])))
df_plot <- rename(df, param=Exp_score)
plot_dist(df_plot, 'Expression', 'graph/CDRH3_distrib_exp.png', expression(bold(Expression~score)))

df_plot <- rename(df, param=motif_CDRH3)
plot_motif_analysis(df_plot, 'graph/CDRH3_motif_analysis_whole.png')
df_plot <- rename(df, param=motif_98)
plot_motif_analysis(df_plot, 'graph/CDRH3_motif_analysis_98.png')
