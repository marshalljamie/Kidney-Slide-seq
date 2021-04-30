import pandas as pd
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
np.random.seed(111)
import argparse

parser = argparse.ArgumentParser()

### Input arguments ###
parser.add_argument('--array_id',help='Unique id for array to be analyzed')
parser.add_argument('--section',help='cortex or medulla')
parser.add_argument('--specimen',help='mouse or human')

args = parser.parse_args()

array_id = args.array_id
section = args.section
specimen = args.specimen

### Variable inits ###
# Set nearest neighbor radius to 25 pixels
radius = 25

# Define sets of cell types to track, specific to specimen of interest
if specimen == 'mouse':
    cell_names = ['PCT_1','PCT_2','EC','MC','Fibroblast','TAL','DCT','CD-IC','CD-PC','GC','MD','Podocyte','Macrophage','Other_Immune','vSMC']
elif specimen == 'human':
    cell_names = ['PCT','EC','MC','Fibroblast','TAL','DCT','CD-IC','CD-PC','GC','MD','Podocyte','Macrophage','Other_Immune','vSMC']

# input_path is path to file with beads x features for all curated cell type calls in an array
# features = {'barcode','x','y','cell_type','section'}
input_path = '{array_id}_allcells.csv'.format(array_id=array_id)
allcells_info = pd.read_csv(input_path,index_col=0)
allcells_info = allcells_info[allcells_info['section']==section]


# input_path is path to file with beads x features for all Trem2-expressing macrophages in an array
# features = {'barcode','x','y','section'}
input_path = '{array_id}_macrophage_Trem2_info.csv'.format(array_id=array_id)
macro_Trem2_info = pd.DataFrame()
if path.exists(input_path):
    macro_Trem2_info = pd.read_csv(input_path,index_col=0)
    macro_Trem2_info = macro_Trem2_info[macro_Trem2_info['section']==section]
    macro_Trem2_barcodes = list(macro_Trem2_info['barcode'])

# if looking at humans, consider Lyve1-expressing macrophages as well
if specimen == 'human':
    # input_path is path to file with beads x features for all Lyve1-expressing macrophages in an array
    # features = {'barcode','x','y','section'}
    input_path = '{array_id}_macrophage_Lyve1_info.csv'.format(array_id=array_id)
    macro_Lyve1_info = pd.DataFrame()
    if path.exists(input_path):
        macro_Lyve1_info = pd.read_csv(input_path,index_col=0)
        macro_Lyve1_info = macro_Lyve1_info[macro_Lyve1_info['section']==section]
        macro_Lyve1_barcodes = list(macro_Lyve1_info['barcode'])

if specimen == 'mouse':
    gene_names = ['Trem2']
elif specimen == 'human':
    gene_names = ['Trem2','Lyve1']

### Trem2/Lyve1-expressing macrophages nn analysis ###
for it,gene_name in enumerate(gene_names):
    if it == 1:
        cell_names = cell_names[:-1]
    cell_names.append(gene_name)

    immune_info = eval('macro_{}_info'.format(gene_name)).copy()
    if not immune_info.empty:
        # label macrophages expressing gene of interest in allcells_info df
        trem2_indices = np.where(allcells_info['barcode'].isin(eval('macro_{}_barcodes'.format(gene_name))))[0]
        allcells_info.loc[trem2_indices,'cell_type'] = gene_name

        # count total instances of every cell type of interest
        celltype_counts = {}
        for name in cell_names:
            celltype_info = allcells_info[allcells_info['cell_type']==name]
            ct = celltype_info.shape[0]
            celltype_counts[name] = ct

        # compute nearest neighbor graph of all beads
        allcells_coords = np.array(allcells_info[['x','y']])
        nn = radius_neighbors_graph(allcells_coords, radius, mode='connectivity',include_self=False)
        nn = nn.toarray()

        # for every bead, count number of instances of every cell type that is a nearest neighbor
        celltype_nn_cts = {}
        for name in cell_names:
            celltype_nn_cts[name] = 0

        all_nn_info_dict = {}
        for row in range(nn.shape[0]):
            nn_row = nn[row,]
            nn_is_true = np.where(nn_row == 1)[0]
            nn_info = allcells_info.iloc[nn_is_true,]
            for name in cell_names:
                celltype_nn_cts[name] = celltype_nn_cts[name]+nn_info[nn_info['cell_type']==name].shape[0]
            all_nn_info_dict[row] = nn_info

        # find nearest neighbor array of gene-expressing macrophages (called Trem2 regardless if analyzing Trem2 or Lyve1)
        trem2_nn = nn[trem2_indices,] 

        # find nearest neighbors of trem2-expressing immune cells (indices of columns of trem2_nn where value = 1)
        trem2_nn_info_dict = {}
        for row in range(trem2_nn.shape[0]):
            trem2_nn_row = trem2_nn[row,]
            trem2_nn_is_true = np.where(trem2_nn_row == 1)[0]
            trem2_nn_info = allcells_info.iloc[trem2_nn_is_true,]
            trem2_nn_info_dict[row] = trem2_nn_info

        # count number of interactions between Trem2-expressing macrophage and other cells
        # normalize with two methods:
        ## interaction(cell type, Trem2-expressing cell) = 
        ## 1. (number of interactions)/sqrt((tot instances of cell type)*(tot instances of trem2-expressing cell))
        ## 2. (number of interactions)/(tot times cell type is nearest neighbor + tot times trem2-expressing cell is nearest neighbor)
        interactions = pd.DataFrame()
        for i in range(len(trem2_nn_info_dict)):
            instance = trem2_nn_info_dict[i]
            for name in cell_names:
                instance_cell = instance[instance['cell_type']==name]
                num_cell = instance_cell.shape[0]
                num_cell_norm1 = num_cell/np.sqrt(celltype_counts[name]*celltype_counts[gene_name])
                if (celltype_nn_cts[name]+celltype_nn_cts[gene_name]) != 0:
                    num_cell_norm2 = num_cell/(celltype_nn_cts[name]+celltype_nn_cts[gene_name])
                else:
                    num_cell_norm2 = 0
                d = {'trem2_id':[i],'celltype':[name],'interaction_norm1':[num_cell_norm1],'interaction_norm2':[num_cell_norm2]}
                d = pd.DataFrame(d)
                interactions = pd.concat([interactions,d])
        interactions = interactions.fillna(0)
        
        # out_path is path to output file
        out_path = '{array_id}_{section}_macro_{gene_name}_interactions.csv'.format(array_id=array_id,section=section,gene_name=gene_name)
        interactions.to_csv(out_path)





