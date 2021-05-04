library(Seurat)
library(optparse)
library(feather)
set.seed(111)

options(future.globals.maxSize= 5242880000)

option_list = list(
  make_option(c('-f', '--ssobjlst'), type='character', default = NULL, 
              help='list containing seurat objects for all arrays of genotypes to be compared'),
  make_option(c('-u', '--arrayids'), type='character', default = NULL,
              help='.txt file containing list of array ids'),
  make_option(c('-g', '--genos'), type='character', default = NULL,
              help='genotypes of arrays'),
  make_option(c('-n', '--num_genes'), type='integer', default = 4000,
              help='number of variable genes to maintain'),
  make_option(c('-o', '--out'), type='character', default = NULL,
              help='path to output directory')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ssobj_lst = get(load(opt$ssobjlst))
array_ids = read.csv(opt$arrayids,header=FALSE)
array_ids[] = lapply(array_ids,as.character)

for (i in 1:length(ssobj_lst)) {
  DefaultAssay(ssobj_lst[[i]]) <- 'SCT'
}

combined_features <- SelectIntegrationFeatures(object.list = ssobj_lst,nfeatures = opt$num_genes)

cts_lst = list()
ids = list()
for(i in 1:length(ssobj_lst)){
  print(i)
  a = t(as.matrix(ssobj_lst[[i]]@assays$Spatial@counts))
  a = a[,combined_features]
  a = as.data.frame(a)
  cts_lst[[i]] = a
  ids[[i]] = rep(array_ids[i,],nrow(a))
}
cts_lst = do.call(rbind,cts_lst)
cts_df = t(as.matrix(cts_lst))
cts_df = as.data.frame(cts_df)

for(i in 1:length(ids)){
  ids[[i]] = as.data.frame(ids[[i]])
  colnames(ids[[i]]) = c('x')
}
ids = do.call(rbind,ids)

barcodes = colnames(cts_df)
barcodes = as.data.frame(barcodes)
colnames(cts_df) = seq(ncol(cts_df))
dat = cts_df[,colSums(cts_df)>0]
rownames(cts_df) = combined_features

combined_ssobj <- CreateSeuratObject(
  counts = dat,
  project = 'SlideSeq',
  assay = 'Spatial'
)

ids = ids[colSums(cts_df)>0,]
barcodes = barcodes[colSums(cts_df)>0,]
combined_ssobj[['arrayids']] = ids
combined_ssobj[['barcodes']] = barcodes
combined_ssobj = SCTransform(combined_ssobj, verbose = FALSE,assay='Spatial',ncells=3000)
Idents(combined_ssobj) <- 'arrayids'

out_path = file.path(opt$out,paste(opt$genos,'_',toString(opt$num_genes),'_combined.RData',sep=''))
save(combined_ssobj,file=out_path)
