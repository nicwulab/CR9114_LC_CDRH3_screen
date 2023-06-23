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

plot_replicate_cor <- function(df, graphname, param){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2)))
  textsize <- 7
  palette <- c(qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  p <- ggplot(df,aes(x=rep1, y=rep2, color=label, size=log10(input_freq))) +
         geom_point(alpha=0.8, pch=16) +
         scale_color_manual('',values=palette,drop=FALSE) +
         scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-3.5, 0),
                                range = c(0.1, 2.5),
                                breaks = c(-3, -2, -1),
                                labels = c(expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1')))) +
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
         labs(x=bquote(bold(paste(.(param),' (replicate 1)'))),y=bquote(bold(paste(.(param),' (replicate 2)'))))
  ggsave(graphname, p, height=2, width=3)
  }

plot_exp_vs_bind <- function(df, graphname){
  textsize <- 7
  palette <- c(qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  p <-  ggplot(df, aes(x=exp_score, y=bind_score, color=label, size=log10(input_freq))) +
          geom_point(alpha=0.8, pch=16) +
          scale_color_manual('',values=palette,drop=FALSE) +
          scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-3.5, 0), 
                                range = c(0.1, 2.5), 
                                breaks = c(-3, -2, -1), 
                                labels = c(expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1')))) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=1,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          xlim(-0.5,1.5) +
          ylim(-0.5,1.5) +
          ylab('Binding score') +
          xlab("Expression score")
  ggsave(graphname,p,width=2.7, height=2, bg='white', dpi=1200)
  }

plot_gene_vs_bind <- function(df, w, h, ymin, ymax, ylab, graphname){
  textsize <- 7
  palette <- c(qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  p <-  ggplot(df, aes(x=gene, y=param, color=label, group=gene, size=log10(input_freq))) +
          geom_boxplot(width=0.7, outlier.shape=NA, size=0.3, color='black') +
          geom_sina(alpha=0.8, pch=16, maxwidth=0.4) +
          scale_color_manual('',values=palette,drop=FALSE) +
          scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-3.5, 0),
                                range = c(0.1, 2.5),
                                breaks = c(-3, -2, -1),
                                labels = c(expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1')))) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylim(ymin,ymax) +
          ylab(ylab) +
          xlab("")
  ggsave(graphname,p,width=w, height=h, bg='white', dpi=1200)
  }

plot_CDRL3len_vs_bind <- function(df, w, h, ymin, ymax, ylab, graphname){
  textsize <- 7
  palette <- c(qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  p <-  ggplot(df, aes(x=gene, y=param, color=label, group=gene, size=log10(input_freq))) +
          geom_boxplot(width=0.7, outlier.shape=NA, size=0.3, color='black') +
          geom_sina(alpha=0.8, pch=16, maxwidth=0.4) +
          scale_color_manual('',values=palette,drop=FALSE) +
          scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-3.5, 0),
                                range = c(0.1, 2.5),
                                breaks = c(-3, -2, -1),
                                labels = c(expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1')))) +
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
          ylim(ymin,ymax) +
          ylab(ylab) +
          xlab("CDR L3 length")
  ggsave(graphname,p,width=w, height=h, bg='white', dpi=1200)
  }

plot_motif_vs_bind <- function(df, w, h, ymin, ymax, ylab, graphname){
  textsize <- 7
  dodge_value <- 0.65
  palette <- brewer.pal(3,"Set1")
  legend_title <- 'CDR L3 residues 91/96'
  p <-  ggplot(df, aes(x=gene, y=param, fill=motif_type, color=motif_type, size=log10(input_freq))) +
          geom_boxplot(width=0.5, alpha=0.4, outlier.shape=NA, size=0.3, color='black', position=position_dodge(dodge_value)) +
          geom_sina(alpha=0.8, pch=16, position=position_dodge(dodge_value), maxwidth=0.4) +
          scale_color_manual(legend_title,values=palette,drop=FALSE) +
          scale_fill_manual(legend_title,values=palette,drop=FALSE) +
          scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-3.5, 0), 
                                range = c(0.1, 2.5), 
                                breaks = c(-3, -2, -1), 
                                labels = c(expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1')))) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.height=unit(0.8,"line"),
                legend.key.size = unit(0.4,"line")) +
          ylim(ymin,ymax) +
          ylab(ylab) +
          xlab("")
  ggsave(graphname,p,width=w, height=h, bg='white', dpi=1200)
  }

plot_SHM_vs_bind <- function(df, w, h, ymin, ymax, ylab, title, graphname){
  textsize <- 7
  palette <- brewer.pal(3,"Set1")
  legend_title <- 'CDR L3 residues 91/96'
  print (graphname)
  for (motif in c('[F/Y/W]/non-[F/Y/W]', 'Others')){
    cor_test <- cor.test(filter(df, motif_type==motif)$SHM_light_count, filter(df, motif_type==motif)$param, method='pearson')
    print (motif)
    print (paste('R =', cor_test$estimate))
    print (paste('P =', cor_test$p.value))
    }
  p <-  ggplot(df, aes(x=SHM_light_count, y=param, group=motif_type, fill=motif_type, color=motif_type, size=log10(input_freq))) +
          geom_point(alpha=0.8, pch=16) +
          geom_smooth(method = "lm", se = TRUE, size=0.5) +
          scale_color_manual(legend_title,values=palette,drop=FALSE) +
          scale_fill_manual(legend_title,values=palette,drop=FALSE) +
          scale_size_continuous('Occurrence\nfrequency',
                                limits = c(-3.5, 0),
                                range = c(0.1, 2.5),
                                breaks = c(-3, -2, -1),
                                labels = c(expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1')))) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=1,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.height=unit(1,"line"),
                legend.key.size = unit(0.4,"line")) +
          ggtitle(title) +
          ylim(ymin,ymax) +
          ylab(ylab) +
          xlab(expression(bold(number~of~V['L']~SHMs)))
  ggsave(graphname,p,width=w, height=h, bg='white', dpi=1200)
  }

labeling_by_antigen <- function(ID, antigen, epi){
  if (grepl('stop',ID)){return ('Nonsense')}
  else if (antigen=='Influenza'){return ('Influenza')}
  else {return ('Others')}
  }

labeling_by_IGLV <- function(ID, Vgene){
  if (grepl('stop',ID)){return ('Nonsense')}
  if (grepl('IGLV1',Vgene)){return ('IGLV1')}
  if (grepl('IGLV2',Vgene)){return ('IGLV2')}
  if (grepl('IGLV3',Vgene)){return ('IGLV3')}
  }

motif_classification <- function(motif){
  motif_type <- c()
  for (aa in str_split(motif,'')[[1]]){
    if (aa %in% c('F','Y','W')){
      motif_type <- c(motif_type, '[F/Y/W]')
      }
    else {
      motif_type <- c(motif_type, 'non-[F/Y/W]')
      }
    }
  motif_type <- paste(motif_type, collapse='/')
  if (motif_type != '[F/Y/W]/non-[F/Y/W]'){motif_type <- 'Others'}
  return (motif_type)
  }

#####MAIN######
set.seed(5)
filename  <- 'result/LC_count_table_score.csv'
graphname1 <- 'graph/LC_replicate_cor_exp.png'
graphname2 <- 'graph/LC_replicate_cor_bind.png'
graphname3 <- 'graph/LC_bind_vs_exp.png'
graphname4 <- 'graph/LC_bind_vs_CDRL3_motif_type.png'
graphname5 <- 'graph/LC_exp_vs_CDRL3_motif_type.png'
graphname6 <- 'graph/LC_bind_vs_CDRL3_length.png'
graphname7 <- 'graph/LC_bind_vs_SHM_count_IGLV1.png'
graphname8 <- 'graph/LC_bind_vs_SHM_count_IGLV2.png'
graphname9 <- 'graph/LC_bind_vs_SHM_count_IGLV3.png'
graphname10 <- 'graph/LC_exp_vs_CDRL3_length.png'

df <- read.csv(filename) %>%
        filter(input_freq > 0.002) %>%
        mutate(label = mapply(labeling_by_antigen, ID, antigen, epi)) %>%
        mutate(label=factor(label, levels=c('Influenza', 'Others', 'Nonsense'))) %>%
        mutate(motif_type=mapply(motif_classification, CDRL3_motif)) %>%
        mutate(IGLV_class = mapply(labeling_by_IGLV, ID, Vgene))
print (mean(filter(df,label=='Influenza')$bind_score))
print (mean(filter(df,label=='Others')$bind_score))
print (t.test(filter(df,label=='Influenza')$bind_score,filter(df,label=='Others')$bind_score))
print (mean(filter(df,label=='Influenza')$exp_score))
print (mean(filter(df,label=='Others')$exp_score))
print (t.test(filter(df,label=='Influenza')$exp_score,filter(df,label=='Others')$exp_score))
print (filter(df,label=='Influenza', bind_score>0.8)$ID)
df_exp  <- df %>%
             rename(rep1=exp_score_1) %>%
             rename(rep2=exp_score_2)
df_bind <- df %>%
             rename(rep1=bind_score_1) %>%
             rename(rep2=bind_score_2)
plot_replicate_cor(df_exp, graphname1, 'Expression score')
plot_replicate_cor(df_bind, graphname2, 'Binding score')
plot_exp_vs_bind(df, graphname3)

df <- df %>% mutate(param=bind_score)
plot_CDRL3len_vs_bind(mutate(filter(df,label!='Nonsense'),gene=CDRL3_len), 3, 2, -0.5, 1.5, 'Binding score', graphname6)
plot_motif_vs_bind(mutate(filter(df,label!='Nonsense'),gene=IGLV_class), 3.5, 2, -0.5, 1.5, 'Binding score', graphname4)
plot_SHM_vs_bind(filter(df,label!='Nonsense' & IGLV_class=='IGLV1'), 3, 2, -0.5, 1.5, 'Binding score', 'IGLV1', graphname7)
plot_SHM_vs_bind(filter(df,label!='Nonsense' & IGLV_class=='IGLV2'), 3, 2, -0.5, 1.5, 'Binding score', 'IGLV2', graphname8)
plot_SHM_vs_bind(filter(df,label!='Nonsense' & IGLV_class=='IGLV3'), 3, 2, -0.5, 1.5, 'Binding score', 'IGLV3', graphname9)

df <- df %>% mutate(param=exp_score)
plot_motif_vs_bind(mutate(filter(df,label!='Nonsense'),gene=IGLV_class), 3, 2, 0.5, 1.5, 'Expression score', graphname5)
plot_CDRL3len_vs_bind(mutate(filter(df,label!='Nonsense'),gene=CDRL3_len), 3, 2, -0.5, 1.5, 'Expression score', graphname10)


#T-test
df <- filter(df, label!='Nonsense') 
print ('IGLV1')
df_test <- filter(df, IGLV_class=='IGLV1')
print (t.test(filter(df_test,motif_type=='[F/Y/W]/non-[F/Y/W]')$bind_score, filter(df_test, motif_type=='Others')$bind_score)$p.value)
print (t.test(filter(df_test,motif_type=='[F/Y/W]/non-[F/Y/W]')$exp_score, filter(df_test, motif_type=='Others')$exp_score)$p.value)
print ('IGLV2')
df_test <- filter(df, IGLV_class=='IGLV2')
print (t.test(filter(df_test,motif_type=='[F/Y/W]/non-[F/Y/W]')$bind_score, filter(df_test, motif_type=='Others')$bind_score)$p.value)
print (t.test(filter(df_test,motif_type=='[F/Y/W]/non-[F/Y/W]')$exp_score, filter(df_test, motif_type=='Others')$exp_score)$p.value)
print ('IGLV3')
df_test <- filter(df, IGLV_class=='IGLV3')
print (t.test(filter(df_test,motif_type=='[F/Y/W]/non-[F/Y/W]')$bind_score, filter(df_test, motif_type=='Others')$bind_score)$p.value)
print (t.test(filter(df_test,motif_type=='[F/Y/W]/non-[F/Y/W]')$exp_score, filter(df_test, motif_type=='Others')$exp_score)$p.value)

print ('CDR L3 length')
print (t.test(filter(df,CDRL3_len==11)$bind_score, filter(df,CDRL3_len==13)$bind_score)$p.value)
print (t.test(filter(df,CDRL3_len==11)$bind_score, filter(df,CDRL3_len==12)$bind_score)$p.value)
print (t.test(filter(df,CDRL3_len==11)$bind_score, filter(df,CDRL3_len==11)$bind_score)$p.value)
print (t.test(filter(df,CDRL3_len==11)$bind_score, filter(df,CDRL3_len==10)$bind_score)$p.value)
print (t.test(filter(df,CDRL3_len==11)$bind_score, filter(df,CDRL3_len==9)$bind_score)$p.value)
