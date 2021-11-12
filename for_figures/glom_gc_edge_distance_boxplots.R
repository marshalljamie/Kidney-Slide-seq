library(ggplot2)

# input_path is path to file with glomeruli x features for glomeruli in all DKD (BTBR ob/ob) and WT (BTBR wt/wt) arrays
# features = {'min_d','geno','array_id'}, where min_d is minimum distance between edge of glom-gc pair, geno is either DKD or WT,
# array_id is unique id of array that glomerulus originates
input_path = 'glom_gc_distance_dat.csv'
dat = read.csv(input_path)
dat['X'] = NULL
dat$pheno <- factor(dat$pheno,
                    levels = c('WT','DKD'),ordered = TRUE)
my_pal = c("#3690C0","#E31A1C")
out_path = 'glom_gc_distance_boxplot.pdf'
cairo_pdf(out_path, height = 7, width = 5,
          family = "Arial", fallback_resolution=300)
print(ggplot(dat, aes(x=pheno, y=min_d, fill=factor(pheno,labels=c('WT','DKD')))) +
        geom_boxplot(outlier.shape=NA,aes(x=factor(pheno,labels=c('WT','DKD')),fill=factor(pheno,labels=c('WT','DKD'))),show.legend = FALSE)+
        geom_jitter(shape=16, position=position_jitter(0.2),show.legend=FALSE)+
        scale_fill_manual(values=my_pal)+
        ylab('Glom-GC distance (pixels)')+
        xlab('')+
        scale_x_discrete(labels= c(expression("BTBR "~italic(wt/wt)),expression("BTBR "~italic(ob/ob))))+
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
