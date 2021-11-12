library(ggplot2)
library(viridis)
library(data.table)
library(optparse)
library(extrafont)
library(stringr)
font_import()
loadfonts()

# input_path is path to file listing unique ids of arrays to be processed
input_path = 'to_process.txt' 
array_ids = read.csv(input_path,header=FALSE)
array_ids[] = lapply(array_ids,as.character)

for(j in 1:nrow(array_ids)){
  array_id = array_ids[j,]
  print(array_id)
  
  # input path is path to file containing beads x features for all beads in array
  # features = {'barcode','x','y'}
  input_path = paste(array_id,'coords.csv',sep='_')
  coords = read.csv(input_path)
  coords['X'] = NULL
  
  # input path is path to file containing beads x features for all beads in array
  # (thresholded cell loading matrix output by NMFreg)
  # features = {'barcode','x','y','max_celltype','cell_type_1',...,'cell_type_n','maxval','thresh_ct'} (loading of each cell type across all beads, max loading per bead, cell type associated with max loading)
  input_path = paste(array_id,'nmf_loadings.csv',sep='_')
  nmf_loadings = read.csv(input_path)
  nmf_loadings['X'] = NULL
  nmf_loadings[['Unnamed..0']] = NULL
  
  # find max celltypes called by nmf
  nmf_max_cell_calls_ind = nmf_loadings[c('max_cell_type')]
  nmf_loadings['max_cell_type'] = NULL
  celltypes = colnames(nmf_loadings)
  celltypes = as.data.frame(celltypes,stringsAsFactors=FALSE)
  celltypes = celltypes[4:(ncol(nmf_loadings)-2),]
  
  nmf_max_cell_calls = list()
  for(i in 1:nrow(nmf_max_cell_calls_ind)){
    nmf_max_cell_calls[[i]] = celltypes[nmf_max_cell_calls_ind[i,]]
  }
  nmf_max_cell_calls = do.call(rbind,nmf_max_cell_calls)
  nmf_loadings[['max_cell_type']] = nmf_max_cell_calls
  
  for(i in 1:length(celltypes)){ 
    print(celltypes[i])
    nmf_onecell_loadings = nmf_loadings[nmf_loadings$max_cell_type == celltypes[i],]
    nmf_onecell_loadings = nmf_onecell_loadings[c('x','y','barcode',celltypes[i])]
    colnames(nmf_onecell_loadings) = c('x','y','barcode','loading')
    nmf_onecell_loadings = nmf_onecell_loadings[c('barcode','x','y','loading')]
    nmf_non_onecell_loadings = subset(coords,!('barcode' %in% nmf_onecell_loadings$barcode))
    nmf_non_onecell_loadings[['loading']] = rep(0.0,nrow(nmf_non_onecell_loadings))
    nmf_final = rbind(nmf_non_onecell_loadings,nmf_onecell_loadings)
    
    out_path=paste(array_id,'nmf',celltypes[i],'maxcalls.pdf',sep='_')
    ggplot(nmf_final, aes(x=x, y=y)) +
      geom_point(size=0.5,aes(colour=loading))+
      scale_colour_viridis(direction=-1)+
      xlab('')+
      ylab('')+
      theme(text = element_text(size=20, colour = 'black',family='Arial'),
            legend.text = element_text(size=20, colour = 'black',family='Arial'),
            legend.title=element_blank(),
            panel.border = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1,
            legend.position = 'none')
    ggsave(out_path, units="in", dpi=300,width=3,height=3)
  }
}
  
for(j in 1:nrow(array_ids)){
  array_id = array_ids[j,]
  print(array_id)
  
  # input path is path to file containing beads x features for all beads in array
  # features = {'barcode','x','y'}
  input_path = paste(array_id,'coords.csv',sep='_')
  coords = read.csv(input_path)
  coords['X'] = NULL
  
  # input_path is path to file with beads x features for all beads in array
  # (loading matrix output by SeuratV3)
  # features = {'cell_type_1',..,'cell_type_n','max'} (loading of each cell type across all beads, and max loading per bead)
  input_path = paste(array_id,'seurat_loadings.csv',sep='_')
  seurat_loadings = read.csv(input_path,row.names=1)
  celltypes = colnames(seurat_loadings)
  celltypes = celltypes[1:(ncol(seurat_loadings)-1)]
  
  max_cell_loadings = seurat_loadings[['max']]
  max_cell_calls = list()
  for(i in 1:nrow(seurat_loadings)){
    col = seurat_loadings[i,]
    max_cell_calls[[i]] = celltypes[min(which(col==max(col)))]
  }
  max_cell_calls = do.call(rbind,max_cell_calls)
  seurat_loadings[['max_cell_type']] = max_cell_calls
  seurat_loadings[['barcode']] = rownames(seurat_loadings)
  
  for(i in 1:length(celltypes)){
    seurat_onecell_loadings = seurat_loadings[seurat_loadings$max_cell_type == celltypes[i],]
    seurat_onecell_loadings = merge(seurat_onecell_loadings,coords,by='barcode')
    seurat_onecell_loadings = seurat_onecell_loadings[c('barcode','x','y',celltypes[i])]
    colnames(seurat_onecell_loadings) = c('barcode','x','y','loading')
    
    seurat_non_onecell_loadings = subset(coords,!('barcode' %in% seurat_onecell_loadings$barcode))
    seurat_non_onecell_loadings[['loading']] = rep(0.0,nrow(seurat_non_onecell_loadings))
    seurat_final = rbind(seurat_non_onecell_loadings,seurat_onecell_loadings) 
    
    out_path=paste(array_id,'_seurat_',celltypes[i],'_maxcalls.pdf',sep='')
    ggplot(seurat_final, aes(x=x, y=y)) +
      geom_point(size=0.5,aes(colour=loading))+
      scale_colour_viridis(direction=-1)+
      xlab('')+
      ylab('')+
      theme(text = element_text(size=20, colour = 'black',family='Arial'),
            legend.text = element_text(size=20, colour = 'black',family='Arial'),
            legend.title=element_blank(),
            panel.border = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1,
            legend.position = 'none')
    ggsave(out_path, units="in", dpi=300,width=3,height=3)
  } 
}