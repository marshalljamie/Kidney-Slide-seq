library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)

# input_path is path to file with disease groups (e.g., healthy vs injured humans) x cell types; each element of the matrix is the average interaction frequency of fixed control 
#cell type (e.g. Macrophages in medulla) and cell type specified by column in disease group 
input_path = 'human_Macrophage_medulla_interactions_avg_by_grp.csv'
dat = read.csv(input_path)
row.names(dat) = dat[['X']]
dat['X'] = NULL
dat = t(dat)
dat = as.data.frame(dat)
dat$celltype = row.names(dat)
dat = dat %>% 
  mutate(celltype = str_replace(celltype, "CD.IC", "CD-IC"))
dat = dat %>% 
  mutate(celltype = str_replace(celltype, "CD.PC", "CD-PC"))

# input_path is path to file with celltype x statistics (pval, pval_adj, significance status)
# where statistics reflect significance based on mann whitney u test comparing all interaction frequencies of fixed control cell type 
# and cell type indicated by row between two disease groups
input_path = 'human_medulla_Macrophage_mannu_pvals.csv'
stats_dat = read.csv(input_path)
stats_dat[c('X','pval','sig_stat')] = NULL

dat = merge(dat,stats_dat,by='celltype')
dat$log_p = -log10(dat$pval_adj)
dat['pval_adj'] = NULL
temp = melt(dat, id=c("celltype","log_p"))
temp_filt = temp[temp$log_p > -log10(0.05),]

out_path = 'human_Macrophage_medulla_interactions_avg_by_grp.pdf'
cairo_pdf(out_path, height = 8, width = 7,
          family = "Arial", fallback_resolution=300)
labels = unique(temp$celltype)
labels = rev(labels)
labels[match('Trem2',labels)] = expression(italic(Trem2)~"+")
labels[match('Lyve1',labels)] = expression(italic(LYVE1)~"+")
print(
  ggplot(temp,aes(x=variable, y = celltype, color = value, size = log_p )) + 
    geom_point(dat =temp_filt ,colour='black',pch=21,fill=NA,stroke=2)+
    geom_point() +
    scale_size_continuous(range=c(3,8))+
    scale_colour_gradient(low='blue',high='red',name='Avg interaction\n frequency')+
    labs(size="-log10(adj pval)")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20, colour = 'black',family='Arial'),
          axis.text = element_text(size=20, colour =  'black',family='Arial'),
          legend.text = element_text(size=20, colour = 'black',family='Arial'),
          legend.key.height = unit(0.8, 'cm'),
          axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5),
          legend.key = element_rect(fill = "white"))+
    scale_y_discrete(labels=labels)+
    xlab('')+
    ylab('')
)
dev.off()