library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(feather)
library(optparse)
library(arrow)
set.seed(111)

option_list = list(
  make_option(c('-i', '--inputDir'), type='character', default = NULL, help = 'Input path'),
  make_option(c('-p', '--arrayid'), type='character', default = NULL, help = 'Unique id of array to be analyzed'),
  make_option(c('-o', '--outDir'), type='character', default = NULL, help = 'Output path')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#### Begin: create Seurat object ####
array_ids = read.csv(opt$arrayids,header=FALSE)
array_ids[] = lapply(array_ids,as.character)

ssobj_list = list()
for (i in 1:nrow(array_ids)){
  array_id = array_ids[i,]
  
  # input path is path to file with beads x features
  # features = {'barcode','x','y'}
  input_path = paste(array_id,'coords.csv',sep='_')
  crds <- read.csv(input_path)
  # input path is path to file with beads x features
  # features = {'barcode','gene1,'gene2',...,'genen'}
  input_path = paste(array_id,'counts.feather',sep='_')
  cts<- arrow::read_feather(input_path)
  
  # check this
  cts['X'] = NULL
  crds['X'] = NULL
  
  barcodes = cts['barcode']
  barcodes = as.vector(barcodes)
  cts['barcode'] = NULL
  cts = as.data.frame(cts)
  rownames(cts) = barcodes[['barcode']]
  cts = t(cts)
  
  ssobj <- CreateSeuratObject(
    counts = cts,
    project = 'SlideSeq',
    assay = 'Spatial'
  )
  
  crds['barcode'] = NULL
  colnames(crds) = c('xcoord','ycoord')
  rownames(crds) = barcodes[['barcode']]
  
  ssobj[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = crds
  )
  #### End: create Seurat object ####
  
  ssobj[['barcode']] = rownames(ssobj@meta.data)
  ssobj[['array_id']] = rep(puck_id,nrow(ssobj@meta.data))
  ssobj[['log_nCount_Spatial']] <- log(ssobj[['nCount_Spatial']])
  ssobj <- NormalizeData(ssobj,normalization.method='LogNormalize',scale.factor = 10000)
  ssobj <- SCTransform(ssobj, assay = "Spatial", ncells = 3000, verbose = FALSE)
  ssobj <- RunPCA(ssobj)
  ssobj <- RunUMAP(ssobj, dims = 1:30)
  ssobj <- FindNeighbors(ssobj, dims = 1:30)
  ssobj <- FindClusters(ssobj, resolution = 0.3, verbose = FALSE)
  ssobj <- FindVariableFeatures(ssobj, selection.method = "vst", nfeatures = 2000)
  
  # out_path is path to output file
  out_path=paste(array_id,'_ssobj.Robj',sep='')
  save(ssobj,file=out_path)
}

