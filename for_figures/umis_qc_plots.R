library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(extrafont)
library(stringr)
library(purrr)
font_import()
loadfonts()

### umis/bead histograms in mice ###
# input_path is path to file with instances of beads x features for all mouse arrays of interest
# features = {'umis','geno'}, where umis are log10-transformed; geno can be WT/UMOD-KI or WT/DKD (a.k.a. BTBR wt/wt/BTBR ob/b)
geno1 = 'UMOD-WT'
geno2 = 'UMOD-KI'
input_path = paste(geno1,geno2,'umis_per_bead.csv',sep='_')
umis_per_bead = read.csv(input_path)
umis_per_bead['X'] = NULL

x = umis_per_bead
x$pheno <- factor(x$geno,
                  levels = c(geno2,geno1),ordered = TRUE)
median <- x %>% 
  pull(umis) %>% 
  median() %>%
  signif(6)

out_path = paste(paste(geno2,geno1,'umis_per_bead_hist',sep='_'),'.pdf',sep='')
cairo_pdf(out_path, height = 5, width = 7,
          family = "Arial", fallback_resolution=300)
print(
  ggplot(x,aes(x=umis, fill=geno)) +
    geom_histogram(position = 'identity',bins=100,alpha=0.5) +
    # include below if comparing WT/DKD
    #scale_fill_manual(values=c("red", "#404080"),labels = c(expression("BTBR "~italic(ob/ob)),expression("BTBR "~italic(wt/wt)))) + 
    scale_fill_manual(values=c("red", "#404080"),labels = c('UMOD_KI','WT')) + 
    scale_colour_manual(values=c("red", "#404080")) +
    xlab('log10(UMIs/bead)')+
    ylab('Frequency')+
    theme_bw()+
    coord_cartesian(xlim = c(0,5),expand=FALSE)+
    theme(text = element_text(size=20, colour = 'black',family='Arial'),
          axis.text = element_text(size=20, colour =  'black',family='Arial'),
          legend.text = element_text(size=20, colour = 'black',family='Arial'),
          legend.title=element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.75, "cm"),
          legend.key.width = unit(0.85, "cm"))+
    guides(fill = guide_legend(reverse=TRUE))+
    geom_vline(xintercept=median,color="black", linetype="dashed", size=0.75)
)
dev.off()

### genes/bead histograms in mice ###
# input_path is path to file with instances of beads x features for all mouse arrays of interest
# features = {'genes','geno'}, where genes are log10-transformed; geno can be WT/UMOD-KI or WT/DKD (a.k.a. BTBR wt/wt/BTBR ob/b)
input_path = paste(geno1,geno2,'genes_per_bead.csv',sep='_')
genes_per_bead = read.csv(input_path)
genes_per_bead['X'] = NULL

x = genes_per_bead
x$pheno <- factor(x$geno,
                  levels = c(geno2,geno1),ordered = TRUE)
median <- x %>% 
  pull(genes) %>% 
  median() %>%
  signif(6)

out_path = paste(paste(geno2,geno1,'genes_per_bead_hist',sep='_'),'.pdf',sep='')
cairo_pdf(out_path, height = 5, width = 7,
          family = "Arial", fallback_resolution=300)
print(
  ggplot(x,aes(x=genes, fill=geno)) +
    geom_histogram(position = 'identity',bins=100,alpha=0.5) +
    # include below if comparing WT/DKD
    #scale_fill_manual(values=c("red", "#404080"),labels = c(expression("BTBR "~italic(ob/ob)),expression("BTBR "~italic(wt/wt)))) +
    scale_fill_manual(values=c("red", "#404080"),labels = c('UMOD_KI','WT')) +
    scale_colour_manual(values=c("red", "#404080")) +
    xlab('log10(genes/bead)')+
    ylab('Frequency')+
    theme_bw()+
    coord_cartesian(xlim = c(0,5),expand=FALSE)+
    theme(text = element_text(size=20, colour = 'black',family='Arial'),
          axis.text = element_text(size=20, colour =  'black',family='Arial'),
          legend.text = element_text(size=20, colour = 'black',family='Arial'),
          legend.title=element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.75, "cm"),
          legend.key.width = unit(0.85, "cm"))+
    guides(fill = guide_legend(reverse=TRUE))+
    geom_vline(xintercept=median,color="black", linetype="dashed", size=0.75) 
)
dev.off()

### human umis/bead histogram ###
# input_path is path to file with instances of beads x features for all human arrays
# features = {'umis'}, where umis are log10-transformed
input_path = 'human_umis_per_bead.csv'
umis_per_bead = read.csv(input_path)
umis_per_bead['X'] = NULL

x = umis_per_bead
median <- x %>% 
  pull(umis) %>% 
  median() %>%
  signif(6)

out_path = 'human_umis_per_bead_hist.pdf'
cairo_pdf(out_path, height = 5, width = 7,
          family = "Arial", fallback_resolution=300)
print(
  ggplot(x,aes(x=umis)) +
    geom_histogram(position = 'identity',bins=100,alpha=0.5,fill='#69b3a2') +
    xlab('log10(UMIs/bead)')+
    ylab('Frequency')+
    theme_bw()+
    coord_cartesian(xlim = c(0,5),expand=FALSE)+
    theme(text = element_text(size=20, colour = 'black',family='Arial'),
          axis.text = element_text(size=20, colour =  'black',family='Arial'),
          legend.text = element_text(size=20, colour = 'black',family='Arial'),
          legend.title=element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    guides(fill = guide_legend(reverse=TRUE))+
    geom_vline(xintercept=median,color="black", linetype="dashed", size=0.75) 
)
dev.off()

### human genes/bead histogram ###
# input_path is path to file with instances of beads x features for all human arrays
# features = {'genes'}, where genes are log10-transformed
input_path = 'human_genes_per_bead.csv'
genes_per_bead = read.csv(input_path)
genes_per_bead['X'] = NULL

x = genes_per_bead
median <- x %>% 
  pull(genes) %>% 
  median() %>%
  signif(6)

out_path = 'human_genes_per_bead_hist.pdf'
cairo_pdf(out_path, height = 5, width = 7,
          family = "Arial", fallback_resolution=300)
print(
  ggplot(x,aes(x=genes, fill=pheno)) +
    geom_histogram(position = 'identity',bins=100,alpha=0.5,fill='#69b3a2') +
    xlab('log10(genes/bead)')+
    ylab('Frequency')+
    theme_bw()+
    coord_cartesian(xlim = c(0,5),expand=FALSE)+
    theme(text = element_text(size=20, colour = 'black',family='Arial'),
          axis.text = element_text(size=20, colour =  'black',family='Arial'),
          legend.text = element_text(size=20, colour = 'black',family='Arial'),
          legend.title=element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    guides(fill = guide_legend(reverse=TRUE))+
    geom_vline(xintercept=median,color="black", linetype="dashed", size=0.75)
)
dev.off()

### For umis per cell type box plots
celltypes = c('EC','TAL','Fibroblast','Macrophage','Other_Immune','MC','GC','Podocyte','MD','CD-IC','CD-PC','DCT','vSMC')
for (i in 1:length(celltypes)){
  celltype = celltypes[i]
  print(celltype)
  # input_path is path to file with beads x features for beads of celltype class in all DKD/WT (same as BTBR ob/ob/BTBR wt/wt) arrays
  # features = {'umis','geno'}
  input_path = paste('DKD_WT_umis_per_',celltype,'_per_bead.csv',sep='')
  DKDWT_df = read.csv(input_path)
  DKDWT_df['X'] = NULL
  # input_path is path to file with beads x features for beads of celltype class in all UMOD-KI/WT arrays
  # features = {'umis','geno'}
  input_path = paste('UMODKI_UMODWT_umis_per_',celltype,'_per_bead.csv',sep='')
  UMODKIWT_df = read.csv(input_path)
  UMODKIWT_df['X'] = NULL
  # input_path is path to file with beads x features for beads of celltype class in all human arrays
  # features = {'umis','geno'}; geno is just 'human'
  input_path = input_path = paste('human_umis_per_',celltype,'_per_bead.csv',sep='')
  human_df = read.csv(input_path)
  human_df['X'] = NULL
  
  all_df = rbind(DKDWT_df,UMODKIWT_df,human_df)
  all_df$umi_tots = log10(all_df$umi_tots)
  all_df$geno <- factor(all_df$geno,
                         levels = c('UMOD-WT','UMOD-KI','WT','DKD','human'),ordered = TRUE)
  
  out_path = paste(celltypes[i],'_umis_per_bead.pdf',sep='')
  cairo_pdf(out_path, height = 5, width = 4,
            family = "Arial", fallback_resolution=300)
  print(
    ggplot(all_df, aes(x=geno, y=umi_tots, fill=factor(geno,labels=c('UMOD-WT','UMOD-KI','WT','DKD','human')))) +
      geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c('UMOD-WT','UMOD-KI','WT','DKD','human')),fill=factor(geno,labels=c('UMOD-WT','UMOD-KI','WT','DKD','human'))),show.legend = FALSE)+
      scale_fill_brewer(palette = "Spectral")+
      scale_x_discrete(labels=c('WT','UMOD_KI',expression("BTBR "~italic(wt/wt)),expression("BTBR "~italic(ob/ob)),'human'))+
      ylab('')+
      xlab('')+
      theme(axis.text.x = element_text(angle = 60, hjust=1),
            text = element_text(size=20, colour = 'black',family='Arial'),
            axis.text = element_text(size=20, colour =  'black',family='Arial'),
            legend.text = element_text(size=20, colour = 'black',family='Arial'),
            legend.title=element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key=element_blank())
  )
  dev.off()
}

### for PCT umis/celltype box plot 
# handled differently because PCT1/2 in mice were combined into PCT class
# input_path is path to file with beads x features for beads of PCT_1 class in all DKD/WT (same as BTBR ob/ob/BTBR wt/wt) arrays
# features = {'umis','geno'}
input_path = paste('DKD_WT_umis_per_','PCT_1','_per_bead.csv',sep='')
DKDWT_df1 = read.csv(input_path)
DKDWT_df1['X'] = NULL
# input_path is path to file with beads x features for beads of PCT_2 class in all DKD/WT (same as BTBR ob/ob/BTBR wt/wt) arrays
# features = {'umis','geno'}
input_path = paste('DKD_WT_umis_per_','PCT_2','_per_bead.csv',sep='')
DKDWT_df2 = read.csv(input_path)
DKDWT_df2['X'] = NULL
# input_path is path to file with beads x features for beads of PCT_1 class in all UMOD-KI/WT arrays
# features = {'umis','geno'}
input_path = paste('UMODKI_UMODWT_umis_per_','PCT_1','_per_bead.csv',sep='')
UMODKIWT_df1 = read.csv(input_path)
UMODKIWT_df1['X'] = NULL
# input_path is path to file with beads x features for beads of PCT_2 class in all UMOD-KI/WT arrays
# features = {'umis','geno'}
input_path = paste('UMODKI_UMODWT_umis_per_','PCT_2','_per_bead.csv',sep='')
UMODKIWT_df2 = read.csv(input_path)
UMODKIWT_df2['X'] = NULL
# input_path is path to file with beads x features for beads of PCT class in all human arrays
# features = {'umis','geno'}
celltype = 'PCT'
input_path = paste('human_umis_per_',celltype,'_per_bead.csv',sep='')
human_df = read.csv(input_path)
human_df['X'] = NULL

all_df = rbind(DKDWT_df1,DKDWT_df2,UMODKIWT_df1,UMODKIWT_df2,human_df)
all_df$umi_tots = log10(all_df$umi_tots)
all_df$geno <- factor(all_df$geno,
                       levels = c('UMOD-WT','UMOD-KI','WT','DKD','human'),ordered = TRUE)
out_path = paste('PCT','_umis_per_bead.pdf',sep='')
cairo_pdf(out_path, height = 5, width = 4,
          family = "Arial", fallback_resolution=300)
print(
  ggplot(all_df, aes(x=geno, y=umi_tots, fill=factor(geno,labels=c('UMOD-WT','UMOD-KI','WT','DKD','human')))) +
    geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c('UMOD-WT','UMOD-KI','WT','DKD','human')),fill=factor(geno,labels=c('UMOD-WT','UMOD-KI','WT','DKD','human'))),show.legend = FALSE)+
    scale_fill_brewer(palette = "Spectral")+
    scale_x_discrete(labels=c('WT','UMOD_KI',expression("BTBR "~italic(wt/wt)),expression("BTBR "~italic(ob/ob)),'human'))+
    ylab('')+
    xlab('')+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          text = element_text(size=20, colour = 'black',family='Arial'),
          axis.text = element_text(size=20, colour =  'black',family='Arial'),
          legend.text = element_text(size=20, colour = 'black',family='Arial'),
          legend.title=element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank())
)
dev.off()