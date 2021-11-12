library(ggplot2)
library(extrafont)
library(RColorBrewer)
font_import()
loadfonts()

### pct in cortex/medulla for humans ###
cell_types = c('EC','TAL','DCT','CD-IC','CD-PC','Fibroblast','vSMC','Macrophage','Other_Immune','PCT')
sections = c('cortex','medulla')

for(i in 1:length(cell_types)){
  for(j in 1:length(sections)){
    cell_type = cell_types[i]
    section = sections[j]
    if(cell_type == 'PCT' & section == 'medulla'){
      next
    }
    # input_path is path to file with instances of arrays x features for all human arrays
    # features = {'geno','pct'}
    input_path = paste('human',cell_type,'pct_in',paste(section,'.csv',sep=''),sep='_')
    dat = read.csv(input_path)
    print(section)
    print(cell_type)
    if(section == 'medulla'){
      dat$geno <- factor(dat$geno,
                          levels = c('1','2','3','4','5','6','DKD','Injured'),ordered = TRUE)
      my_pal = c("#023858","#045A8D","#0570B0","#3690C0","#74A9CF","#A6BDDB",'#E31A1C','#FC4E2A')
      out_path = paste('human',cell_type,'pct_in',paste(section,'.pdf',sep=''),sep='_')
      cairo_pdf(out_path, height = 5, width = 4,
                family = "Arial", fallback_resolution=300)
      print(ggplot(dat, aes(x=geno, y=pct, fill=factor(geno,labels=c('1','2','3','4','5','6','DKD','Injured')))) +
        geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c('1','2','3','4','5','6','DKD','Injured')),fill=factor(geno,labels=c('1','2','3','4','5','6','DKD','Injured'))),show.legend = FALSE)+
        geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5,show.legend=FALSE,colour='black',fill='black')+
        scale_fill_manual(values=my_pal)+
        ylab('')+
        xlab('')+
        scale_x_discrete(labels= c('Healthy 1','Healthy 2','Healthy 3','Healthy 4','Healthy 5','Healthy 6','DKD','Injured'))+
        theme(axis.text.x = element_text(angle = 60, hjust=1),
              text = element_text(size=20, colour = 'black', family = 'Arial'),
              axis.text = element_text(size=20, colour =  'black', family = 'Arial'),
              legend.text = element_text(size=20, colour = 'black', family = 'Arial'),
              legend.title=element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(),
              panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key=element_blank()))
      dev.off()
    }else{
      dat$geno <- factor(dat$geno,
                          levels = c('1','2','3','4','5','6','7','DKD','Injured'),ordered = TRUE)
      my_pal = c("#023858","#045A8D","#0570B0","#3690C0","#74A9CF","#A6BDDB",'#D0D1E6','#E31A1C','#FC4E2A')
      out_path = paste('human',cell_type,'pct_in',paste(section,'.pdf',sep=''),sep='_'))
      cairo_pdf(out_path, height = 5, width = 4,
                family = "Arial", fallback_resolution=300)
      print(ggplot(dat, aes(x=geno, y=pct, fill=factor(geno,labels=c('1','2','3','4','5','6','7','DKD','Injured')))) +
        geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c('1','2','3','4','5','6','7','DKD','Injured')),fill=factor(geno,labels=c('1','2','3','4','5','6','7','DKD','Injured'))),show.legend = FALSE)+
        geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5,show.legend=FALSE,colour='black',fill='black')+
        scale_fill_manual(values=my_pal)+
        ylab('')+
        xlab('')+
        scale_x_discrete(labels= c('Healthy 1','Healthy 2','Healthy 3','Healthy 4','Healthy 5','Healthy 6','Healthy 7','DKD','Injured'))+
        theme(axis.text.x = element_text(angle = 60, hjust=1),
              text = element_text(size=20, colour = 'black', family = 'Arial'),
              axis.text = element_text(size=20, colour =  'black', family = 'Arial'),
              legend.text = element_text(size=20, colour = 'black', family = 'Arial'),
              legend.title=element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(),
              panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key=element_blank()))
      dev.off()
    }
  }
}

### pct in gloms for humans ###
cell_types = c('podo','mesangial','endothelial')
for(i in 1:length(cell_types)){
  cell_type = cell_types[i]
  # input_path is path to file with instances of glomeruli x features for all human arrays
  # features = {'geno','proportion_of_glom'}
  input_path = paste(paste('prop',cell_type,'in_glom_human',sep='_'),'.csv',sep='')
  dat = read.csv(input_path)
  dat$geno <- factor(dat$geno,
                      levels = c('1','2','3','4','5','6','7','DKD','Injured'),ordered = TRUE)
  br_pal1 <- brewer.pal(9,"YlOrRd")
  my_pal = c("#023858","#045A8D","#0570B0","#3690C0","#74A9CF","#A6BDDB",'#D0D1E6','#E31A1C','#FC4E2A')
  out_path = paste(paste('human',cell_type,'pct_in_glom',sep='_'),'.pdf',sep='')
  cairo_pdf(out_path, height = 5, width = 4,
            family = "Arial", fallback_resolution=300)
  print(ggplot(dat, aes(x=geno, y=pct, fill=factor(geno,labels=c('1','2','3','4','5','6','7','DKD','Injured')))) +
    geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c('1','2','3','4','5','6','7','DKD','Injured')),fill=factor(geno,labels=c('1','2','3','4','5','6','7','DKD','Injured'))),show.legend = FALSE)+
    geom_jitter(shape=16, position=position_jitter(0.2),show.legend=FALSE)+
    scale_fill_manual(values=my_pal)+
    ylab('')+
    xlab('')+
    scale_x_discrete(labels= c('Healthy 1','Healthy 2','Healthy 3','Healthy 4','Healthy 5','Healthy 6','Healthy 7','DKD','Injured'))+
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          text = element_text(size=20, colour = 'black', family = 'Arial'),
          axis.text = element_text(size=20, colour =  'black', family = 'Arial'),
          legend.text = element_text(size=20, colour = 'black', family = 'Arial'),
          legend.title=element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank()))
  dev.off()
}

### pct in all sections for mice ###
# DKD = BTBR ob/ob, WT = BTBR wt/wt
# geno1/geno2 can be WT/DKD or UMOD-WT/UMOD-KI
geno1 = 'WT'
geno2 = 'DKD'
cell_types = c('EC','TAL','DCT','CD-IC','CD-PC','Fibroblast','vSMC','Macrophage','Other_Immune','PCT_1','PCT_2')
sections = c('cortex','medulla','all')
for(i in 1:length(cell_types)){
  for(j in 1:length(sections)){
    cell_type = cell_types[i]
    section = sections[j]
    if((cell_type == 'PCT_1' | cell_type == 'PCT_2') & (section == 'medulla' or section == 'all')){
      next
    }
    # input_path is path to file with instances of arrays x features for all mouse arrays
    # features = {'geno','pct'}
    input_path = paste(paste('mouse',cell_type,'pct_in',section,'df',sep='_'),'.csv',sep='')
    dat = read.csv(input_path)
    dat = dat[dat$geno %in% c(geno1,geno2),]
    colnames(dat) = c('geno','pct')
    
    dat$geno <- factor(dat$geno,
                        levels = c(geno1,geno2),ordered = TRUE)
    my_pal = c("#3690C0","#E31A1C")
    out_path = paste(paste(geno1,geno2,cell_type,'pct_in',section,sep='_'),'.pdf',sep='')
    cairo_pdf(out_path, height = 5, width = 3,
              family = "Arial", fallback_resolution=300)
    print(ggplot(dat, aes(x=geno, y=pct, fill=factor(geno,labels=c(geno1,geno2)))) +
            geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c(geno1,geno2)),fill=factor(geno,labels=c(geno1,geno2))),show.legend = FALSE)+
            geom_jitter(shape=16, position=position_jitter(0.2),show.legend=FALSE)+
            scale_fill_manual(values=my_pal)+
            ylab('')+
            xlab('')+
            scale_x_discrete(labels= c(expression("BTBR "~italic(wt/wt)),expression("BTBR "~italic(ob/ob))))+ # will need to be changed if comparing WT/UMOD-KI
            theme(axis.text.x = element_text(angle = 60, hjust=1),
                  text = element_text(size=20, colour = 'black', family = 'Arial'),
                  axis.text = element_text(size=20, colour =  'black', family = 'Arial'),
                  legend.text = element_text(size=20, colour = 'black', family = 'Arial'),
                  legend.title=element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  aspect.ratio=2/1,
                  legend.key=element_blank()))
    dev.off()
  }
}

### pct in glom for mice ###
cell_types = c('podo','endothelial','mesangial')
for(i in 1:length(cell_types)){
  cell_type = cell_types[i]
  # input_path is path to file with instances of glomeruli x features for all mouse arrays
  # features = {'geno','proportion_of_glom'}
  input_path = paste(paste('prop',cell_type,'in_glom_mice',sep='_'),'.csv',sep='')
  dat = read.csv(input_path)
  dat = dat[dat$geno %in% c(geno1,geno2),]
  colnames(dat) = c('geno','pct')
  dat$geno <- factor(dat$geno,
                      levels = c(geno1,geno2),ordered = TRUE)
  my_pal = c("#3690C0","#E31A1C")
  out_path = paste(paste(geno1,geno2,cell_type,'pct_in_glom',sep='_'),'.pdf',sep=''))
  cairo_pdf(out_path, height = 5, width = 3,
            family = "Arial", fallback_resolution=300)
  print(ggplot(dat, aes(x=geno, y=pct, fill=factor(geno,labels=c(geno1,geno2)))) +
          geom_boxplot(outlier.shape=NA,aes(x=factor(geno,labels=c(geno1,geno2)),fill=factor(geno,labels=c(geno1,geno2))),show.legend = FALSE)+
          geom_jitter(shape=16, position=position_jitter(0.2),show.legend=FALSE)+
          scale_fill_manual(values=my_pal)+
          ylab('')+
          xlab('')+
          scale_x_discrete(labels= c(expression("BTBR "~italic(wt/wt)),expression("BTBR "~italic(ob/ob))))+ # will need to be changed if comparing WT/UMOD-KI
          theme(axis.text.x = element_text(angle = 60, hjust=1),
                text = element_text(size=20, colour = 'black', family = 'Arial'),
                axis.text = element_text(size=20, colour =  'black', family = 'Arial'),
                legend.text = element_text(size=20, colour = 'black', family = 'Arial'),
                legend.title=element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(),
                panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=2/1,
                legend.key=element_blank()))
  dev.off()
}













