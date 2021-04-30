library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(feather)
library(optparse)
set.seed(111)

option_list = list(
  make_option(c('-p', '--arrayids'), type='character', default = NULL, help = 'Unique ids of arrays to be analyzed'),
  make_option(c('-r', '--reference'), type='character', default = NULL, help = 'Single cell reference object'),
  make_option(c('-o', '--outDir'), type='character', default = NULL, help = 'Output path')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

array_ids = read.csv(opt$arrayids,header=FALSE)
array_ids[] = lapply(array_ids,as.character)

for (i in 1:nrow(array_ids)){
  array_id = array_ids[i,]
  
  ref = get(load(opt$reference))
  # input_path is path to seurat objects for all arrays, after seurat-preprocessing.R
  input_path=paste(array_id,'ssobj.Robj',sep='_')
  ssobj = get(load(input_path))
  
  ssobj[['array']] = rep(array_id,nrow(ssobj@meta.data))
  
  anchors <- FindTransferAnchors(reference = ref, query = ssobj, normalization.method = "SCT", 
                                 npcs = 50)
  predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE, 
                                    weight.reduction = ssobj[["pca"]], dims = 1:50)
  ssobj[["predictions"]] <- predictions.assay
  DefaultAssay(ssobj) <- "predictions"
  
  celltypes = as.data.frame(unique(ref@meta.data['celltype']))
  celltypes = lapply(celltypes[['celltype']],as.character)
  celltypes = do.call(rbind,celltypes)
  celltypes_temp = gsub('_','-',celltypes[,1])
  celltypes_name_adj = gsub('/','-',celltypes[,1])
  
  # generate spatial plots of prediction loadings per cell type
  for(i in 1:length(celltypes)){
    SpatialFeaturePlot(ssobj, features = c(celltypes_temp[i]), alpha = c(0.1, 1)) 
    out_path=file.path(opt$outDir,paste(array_id,'seurat',celltypes_name_adj[i],'loadings.pdf',sep='_'))
    ggsave(out_path, units="in", dpi=300)
  }
  
  # output cell type predictions (cell type corresponding to max prediction loading per bead)
  cell_loadings = ssobj@assays[['predictions']]@data
  max_cell_loadings = cell_loadings[nrow(cell_loadings),]
  max_cell_calls = list()
  for(i in 1:ncol(cell_loadings)){
    col = cell_loadings[,i]
    max_cell_calls[[i]] = celltypes[min(which(col==max(col)))]
  }
  max_cell_calls = do.call(rbind,max_cell_calls)
  
  ssobj[['max_celltype']] = max_cell_calls
  
  out_path=file.path(opt$outDir,paste(array_id,'_ssobj.Robj',sep='_'))
  save(ssobj,file=out_path)
}






