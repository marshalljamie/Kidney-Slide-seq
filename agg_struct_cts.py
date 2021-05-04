import pandas as pd
import numpy as np
import os
from os import path
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--specimen',help='mouse or human')
parser.add_argument('--array_id',help='unique id for array to be analyzed')
parser.add_argument('--section',help='cortex or medulla')

args = parser.parse_args()

specimen = args.specimen
array_id = args.array_id
section = args.section

def aggregate_struct(celltype_counts):    
    clust_to_agg = pd.DataFrame()
    
    name = 'cluster'
    for i in np.unique(np.array(celltype_cts[name])):  
        clust = celltype_counts[celltype_cts[name] == i]
        clust_genesOnly = clust.iloc[:,6:]
        clust_genes_avg = clust_genesOnly.mean(axis=0)
        clust_genes_avg = pd.DataFrame(clust_genes_avg)
        clust_genes_avg = clust_genes_avg.T
        clust_genes_avg['cluster'] = i
        clust_to_agg = pd.concat([clust_to_agg,clust_genes_avg])

# input path is path to file with beads x features for curated and partitioned cell type of interest in section of interest
# features = {'barcode','x','y','cell_type','section','cluster','gene1','gene2',...,'genen'}
input_path = '{array_id}_{cell_type}_{section}_cts.pkl'.format(array_id=array_id,section=section,cell_type=cell_type)
if path.exists(input_path):
    celltype_counts = pd.read_pickle(input_path)
    celltype_agg = aggregate_struct(celltype_counts)

    # out_path is path to output file
    out_path = '{array_id}_{cell_type}_{section}_agg_cts.csv'.format(array_id=array_id,cell_type=cell_type,section=section)
    celltype_agg.to_csv(out_path)