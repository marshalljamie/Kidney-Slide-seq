library(Seurat)
library(optparse)
set.seed(111)

options(future.globals.maxSize= 2097152000)

option_list = list(
  make_option(c('-f', '--ssobjlst'), type='character', default = NULL, 
              help='list containing seurat objects for all arrays of one genotype'),
  make_option(c('-u', '--arrayids'), type='character', default = NULL,
              help='.txt file containing list of array ids'),
  make_option(c('-g', '--geno'), type='character', default = NULL,
              help='genotype of arrays'),
  make_option(c('-o', '--out'), type='character', default = NULL,
              help='path to output directory')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# dictionary of mappings betwee array ids and batch ids
array_to_batch = list(
  '191109_06'= 1,
  '191109_07'= 1,
  '191109_08'= 1,
  '191109_09'= 1,
  '191109_10'= 1,
  '191109_11'= 1,
  '191109_12'= 1,
  '191109_14'= 1,
  '191109_18'= 1,
  '191109_19'= 1,
  '191109_20'= 1,
  '191109_23'= 1,
  '191109_24'= 1,
  '191112_04'= 2,
  '191112_02'= 2,
  '191112_05'= 2,
  '191112_07'= 2,
  '191112_12'= 2,
  '191112_13'= 2,
  '191112_08'= 2,
  '191204_03'= 3,
  '191204_04'= 3,
  '191204_05'= 3,
  '191204_07'= 3,
  '191204_08'= 3,
  '191204_09'= 3,
  '191204_12'= 3,
  '191204_13'= 4,
  '191204_14'= 4,
  '191204_15'= 4,
  '191204_16'= 4,
  '191204_17'= 4,
  '191204_18'= 4,
  '191204_20'= 4,
  '191204_22'= 4,
  '191204_23'= 4,
  '191204_24'= 4,
  '191206_01'= 4,
  '191206_02'= 4,
  '191206_03'= 4,
  '191206_04'= 4,
  '191223_01'= 5,
  '191223_02'= 5,
  '191223_03'= 5,
  '191223_04'= 5,
  '191223_05'= 5,
  '191223_07'= 5,
  '191223_08'= 5,
  '191223_09'= 5,
  '191223_10'= 5,
  '191223_11'= 5,
  '191223_12'= 5,
  '191223_13'= 5,
  '191223_14'= 5,
  '191223_15'= 6,
  '191223_17'= 6,
  '191223_18'= 6,
  '191223_19'= 6,
  '191223_20'= 6,
  '191223_21'= 6,
  '191223_22'= 6,
  '191223_23'= 6,
  '191223_24'= 6,
  '200102_09'= 6,
  '200104_01'= 7,
  '200104_02'= 7,
  '200104_03'= 7,
  '200104_04'= 7,
  '200104_05'= 7,
  '200104_06'= 7,
  '200104_07'= 7,
  '200104_09'= 7,
  '200104_10'= 7,
  '200104_14'= 7,
  '200104_15'= 8,
  '200104_16'= 8,
  '200104_17'= 8,
  '200104_18'= 8,
  '200104_19'= 8,
  '200104_20'= 8,
  '200104_21'= 8,
  '200104_23'= 8,
  '200115_01'= 9,
  '200115_02'= 9,
  '200115_03'= 9,
  '200115_04'= 9,
  '200115_05'= 9,
  '200115_06'= 9,
  '200115_07'= 9,
  '200115_09'= 9,
  '200115_10'= 10,
  '200115_11'= 10,
  '200115_12'= 10,
  '200115_13'= 10,
  '200115_14'= 10,
  '200115_15'= 10,
  '200115_16'= 10,
  '200115_17'= 10,
  '200115_18'= 10,
  '200113_07'= 11,
  '200113_08'= 11,
  '200113_09'= 11,
  '200113_10'= 11,
  '200113_11'= 11,
  '200113_12'= 11,
  '200121_01'= 11,
  '200121_03'= 11,
  '200127_01'= 12,
  '200127_02'= 12,
  '200127_03'= 12,
  '200127_04'= 12,
  '200127_05'= 12,
  '200127_06'= 12,
  '200127_07'= 12,
  '200127_08'= 12,
  '200127_09'= 12,
  '200127_10'= 12,
  '200131_24'= 13,
  '200131_25'= 13,
  '200131_26'= 13,
  '200205_13'= 13,
  '200205_16'= 13,
  '200205_17'= 13,
  '200205_18'= 13,
  '200205_19'= 13,
  '200131_13'= 14,
  '200131_20'= 14,
  '200131_21'= 14,
  '200131_22'= 14,
  '200131_23'= 14,
  '200205_20'= 16,
  '200205_21'= 16,
  '200205_22'= 16,
  '200205_23'= 16,
  '200205_24'= 16,
  '200205_25'= 16,
  '200210_01'= 16,
  '200210_02'= 16,
  '200210_03'= 16,
  '200210_04'= 16
)

geno1_ssobj_lst = get(load(opt$ssobjlstone))
ssobj_lst = get(load(opt$ssobjlst))

array_ids = read.csv(opt$arrayids,header=FALSE)
array_ids[] = lapply(array_ids,as.character)

for (i in 1:length(ssobj_lst)) {
  DefaultAssay(ssobj_lst[[i]]) <- 'SCT'
}

combined_features <- SelectIntegrationFeatures(object.list = combined_ssobj_lst,nfeatures = 2000)

cts_lst = list()
batches = list()
for(i in 1:length(ssobj_lst)){
  print(i)
  a = t(as.matrix(ssobj_lst[[i]]@assays$Spatial@counts))
  a = a[,combined_features]
  a = as.data.frame(a)
  cts_lst[[i]] = a
  batches[[i]] = rep(array_to_batch[[array_ids[i,]]],nrow(a))
}
cts_lst = do.call(rbind,cts_lst)
cts_df = t(as.matrix(cts_lst))
cts_df = as.data.frame(cts_df)

for(i in 1:length(batches)){
  batches[[i]] = as.data.frame(batches[[i]])
  colnames(batches[[i]]) = c('x')
}
batches = do.call(rbind,batches)

colnames(cts_df) = seq(ncol(cts_df))
cts_df = cts_df[,colSums(cts_df)>0]
rownames(cts_df) = combined_features

combined_ssobj <- CreateSeuratObject(
  counts = cts_df,
  project = 'SlideSeq',
  assay = 'Spatial'
)

batches = batches[colSums(cts_df)>0,]
combined_ssobj[['batch']] = batches
combined_ssobj = SCTransform(combined_ssobj, verbose = FALSE,assay='Spatial',ncells=3000)
Idents(combined_ssobj) <- 'batch'

combined_batch_markers = FindAllMarkers(combined_ssobj, assay = 'SCT',test.use = 'negbinom',only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
out_path = file.path(opt$out,paste(opt$geno,'_batch_markers.csv',sep='')) 
combined_batch_markers.to_csv(out_path)









