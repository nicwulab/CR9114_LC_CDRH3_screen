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

plot_aa_freq <- function(df, w, xlab, graphname){
  aa_levels <- sort(unique(df$aa))
  df <- mutate(df, aa=factor(aa, levels=aa_levels))
  colorscale <- brewer.pal(3,"Set1")
  textsize <- 7
  p <- ggplot() +
	  geom_bar(data=df, aes(x=aa,y=freq,fill=bind), stat="identity", width=0.8, position=position_dodge(width=0.8), alpha=0.8) +
	  scale_fill_manual(values=colorscale) +
	  theme_cowplot(12) +
	  theme(axis.title=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"),
		axis.text.x=element_text(angle = 0, hjust = 0.5, size=textsize, vjust=1, face="bold"),
		legend.key.size=unit(0.12,'in'),
		legend.title=element_blank(),
		legend.text=element_text(size=textsize,face="bold"),
		legend.position='top') +
	  #scale_y_continuous(breaks=c(0,2,4,6,8),labels=c('0','2','4','6','8')) +
	  labs(y=expression(bold('Frequency (%)')),x=xlab)
  ggsave(graphname,p,height=1.5,width=w,bg='white')
  }

df <- read_tsv('result/bind_aa_freq.tsv')
plot_aa_freq(filter(df, resi==91), 1.4, 'residue 91', 'graph/bind_aa_freq_91.png')
plot_aa_freq(filter(df, resi==96), 2, 'residue 96', 'graph/bind_aa_freq_96.png')
