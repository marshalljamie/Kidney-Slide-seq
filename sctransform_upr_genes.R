library(Seurat)
library(arrow)
library(SeuratData)
options(future.globals.maxSize= 5242880000)

option_list = list(
  make_option(c('-i', '--array_ids'), type='character', default = NULL,help = '.txt file with list of UMOD-KI and UMOD-WT unique ids'),
  make_option(c('-f','--inFile'), type='character', default = NULL,help = 'data matrix with UPR gene raw counts for all UMOD-KI and UMOD-WT arrays'),
  make_option(c('-d', '--outDir'), type='character', default = NULL)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

array_ids = opt$array_ids
inFile = opt$inFile
outDir = opt$outDir

# More info on inFile:
# data matrix with beads x features for all UMOD-KI and UMOD-WT arrays
# features ={'barcode','genotype','array_id','gene1','gene2',...,'genen'}
# Raw counts are stored per gene, and genes are UPR pathway genes
dat = arrow::read_feather(inFile)
dat['index'] = NULL
info = dat[c('barcode','genotype','array_id')]
dat['barcode'] = NULL
dat['genotype'] = NULL
dat['array_id'] = NULL
dat = as.data.frame(dat)
col_names <- colnames(dat)
dat[col_names] <- sapply(dat[col_names],as.numeric)
dat = t(dat)
colnames(dat) = seq(ncol(dat))

ssobj <- CreateSeuratObject(
  counts = dat,
  project = 'SlideSeq',
  assay = 'Spatial'
)

ssobj[['barcode']] = info[c('barcode')]
ssobj[['array_id']] = info[c('array_id')]
ssobj[['genotype']] = info[c('genotype')]

ssobj <- SCTransform(ssobj, assay = "Spatial",variable.features.n=NULL)

array_ids = read.csv(array_ids,header=FALSE)
array_ids[] = lapply(array_ids,as.character)

Idents(ssobj) <- 'array_id'

for (i in 1:nrow(array_ids)){
  array_id = array_ids[i,]
  array_ssobj = subset(x=ssobj,idents=array_id)
  barcodes = array_ssobj@meta.data['barcode']
  barcodes$barcodes = as.character(barcodes$barcodes)
  barcodes = barcodes$barcodes
  barcodes = gsub('[[:digit:]]+', '', barcodes)
  
  array_sct_cts = as.data.frame(GetAssayData(object = array_ssobj, assay = 'SCT',slot='data'))
  colnames(array_sct_cts) = barcodes
  array_sct_cts['genes'] = rownames(array_sct_cts)
  # out_path is path to output file
  out_path = file.path(outDir,paste(array_id,'UPR_pathway_genes_sct_cts.feather',sep=''))
  write_feather(array_sct_cts,out_path) 
}




