library(Seurat)
library(ggplot2)
library(SeuratData)

# input_path is path to file containing cells x marker genes, where cells are those confidently called a cell type for all human arrays
input_path = 'human_all_gene_markers.csv'
data = read.csv(input_path)
data['X'] = NULL
celltypes = data$cell_type
data['cell_type'] = NULL
data = t(data)
data = as.data.frame(data)

specimen = 'human' # mouse or human

ssobj = CreateSeuratObject(counts=data)
ssobj[['celltype']] = celltypes
Idents(object=ssobj) <- "celltype"

if (specimen == 'mouse'){
  ssobj@active.ident <- factor(ssobj@active.ident,
                              levels=c("Podocyte",
                                       "MC",
                                       "EC",
                                       "GC",
                                       "vSMC",
                                       "Fibroblast",
                                       "Macrophage",
                                       "Other_Immune",
                                       "TAL",
                                       "MD",
                                       "PCT_2",
                                       "PCT_1",
                                       "DCT",
                                       "CD-IC",
                                       "CD-PC"))
}else{
  ssobj@active.ident <- factor(ssobj@active.ident,
                               levels=c("Podocyte",
                                        "EC",
                                        "vSMC",
                                        "GC",
                                        "MC",
                                        "Fibroblast",
                                        "Macrophage",
                                        "Other_Immune",
                                        "TAL",
                                        "MD",
                                        "DCT",
                                        "CD-IC",
                                        "CD-PC",
                                        "PCT"
                               ))
}

DotPlot(ssobj,features=rownames(data),cols=c('blue','red'))+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
        text = element_text(size=14, colour = 'black',family='Arial'),
        axis.text = element_text(size=14, colour =  'black',family='Arial'),
        legend.text = element_text(size=14, colour = 'black',family='Arial')
        )+
  xlab('Marker Gene')+
  ylab('Cell type')
        
        


