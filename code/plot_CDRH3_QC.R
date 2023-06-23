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

plot_replicate_cor <- function(df, graphname, xlab, ylab, axis_min, axis_max){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2, use='pairwise.complete.obs')))
  textsize <- 7
  class_levels <- c('CR9114 mutant', 'Nonsense', 'Influenza (HA stem)', 'Influenza (non-HA stem)', 'Non-influenza', 'WT')
  palette <- c(qualpal(n = 8, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  palette <- c(palette[1:3], palette[7], 'gray', palette[8])
  df <- filter(df, class != 'WT') %>%
          mutate(class=factor(class, levels=class_levels))
  p <- ggplot(df,aes(x=rep1, y=rep2, color=class, size=log10(avg_input_freq))) +
         geom_point(alpha=0.5, pch=16) +
         scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-5.5, -2.5),
                                range = c(0.1, 0.8),
                                breaks = c(-5, -4, -3),
                                labels = c(expression(bold('10'^'-5')), expression(bold('10'^'-4')), expression(bold('10'^'-3')))) +
         scale_color_manual('',values=palette,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         xlim(axis_min, axis_max) +
         ylim(axis_min, axis_max) +
         labs(x=xlab,y=ylab)
  ggsave(graphname, p, height=1.85, width=3.1,dpi=600)
  }

assign_epi <- function(class, epi){
  if (class=='nonsense'){return ('Nonsense')}
  else if (class=='CR9114 mutant'){return ('CR9114 mutant')}
  else if (class=='WT'){return ('WT')}
  else if (epi=='HA stem'){return ('Influenza (HA stem)')}
  else if (class=='influenza'){return ('Influenza (non-HA stem)')}
  else {return ('Non-influenza')}
  }

#####MAIN######
set.seed(5)
graphname1 <- 'graph/CDRH3_replicate_cor_exp.png'
graphname2 <- 'graph/CDRH3_replicate_cor_KD.png'
df      <- read_csv('result/CDRH3_express_table_summary.csv')
df_exp  <- read_csv('result/CDRH3_KD_table_summary_rep.csv')
df_bind <- read_csv('result/CDRH3_KD_table_summary_class.csv')
df <- inner_join(df,df_bind) %>%
  inner_join(df_exp) %>%
  mutate(log10_Kd=log10(Kd)) %>%
  filter((log10_Kd < -8 & p.value < 0.2) | (log10_Kd >= -8)) %>%
  mutate(class = mapply(assign_epi, class, epi)) %>%
  filter(avg_input_freq > 0)
print (log10(range(df$avg_input_freq)))

df_exp <- df %>%
             rename(rep1=Exp_score_1) %>%
             rename(rep2=Exp_score_2)
plot_replicate_cor(df_exp, graphname1,
                    bquote(bold("Expression score (replicate 1)")), bquote(bold("Expression score (replicate 2)")), -0.77, 3.75)
df_KD  <- df %>%
             rename(rep1=Kd.x) %>%
             rename(rep2=Kd.y) %>%
             mutate(rep1=log10(rep1)) %>%
             mutate(rep2=log10(rep2)) 
plot_replicate_cor(df_KD, graphname2,
                   expression(bold(log['10']~K['D']~'(replicate 1)')), expression(bold(log['10']~K['D']~'(replicate 2)')), -10.5, -4.5)
