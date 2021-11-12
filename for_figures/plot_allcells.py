import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from os import path
import argparse

np.random.seed(111)

parser = argparse.ArgumentParser()

parser.add_argument('--array_id',help='Unique id for array to be processed')
parser.add_argument('--specimen',help='mouse or human')

args = parser.parse_args()

specimen = args.specimen
array_id = args.array_id

# input_path is path to file with beads x features for all curated cell type calls in an array
# features = {'barcode','x','y','cell_type','section'}
input_path = '{array_id}_allcells.csv'.format(array_id=array_id)
allcells_df = pd.read_csv(input_path,index_col=0)
allcells_df = allcells_df.reset_index()
allcells_df = allcells_df.drop(columns={'index'})

if specimen == 'mouse':
    allcells_df['cell_type'] = pd.Categorical(allcells_df['cell_type'], ['EC','PCT_1','PCT_2','Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD'])
elif specimen == 'human':
    allcells_df['cell_type'] = pd.Categorical(allcells_df['cell_type'], ['EC','PCT','Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD'])
allcells_df = allcells_df.sort_values("cell_type")

out_path = '{array_id}_allcells.pdf'.format(array_id=array_id))
subplot_kw = dict(xlim=(0, 6000), ylim=(0, 6000), autoscale_on=False)
fig, ax = plt.subplots(subplot_kw=subplot_kw,figsize=(10,10))

if specimen == 'mouse':
    colors = ['lightblue',"plum","pink",'grey','mediumseagreen','gold','wheat','yellowgreen','olivedrab','rebeccapurple','teal','midnightblue','cyan','orangered','magenta']
    texts = ['EC',"PCT_1","PCT_2",'Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD']
    color_d = {
        'Podocyte' : 'midnightblue',
        'MC' : 'cyan',
        'EC' : 'lightblue',
        'PCT_1' : 'plum',
        'PCT_2' : 'pink',
        'CD-IC' : 'olivedrab',
        'CD-PC' : 'yellowgreen',
        'Other_Immune' : 'grey',
        'Macrophage': 'mediumseagreen',
        'GC' : 'orangered',
        'MD' : 'magenta',
        'TAL' : 'gold',
        'DCT' : 'wheat',
        'vSMC': 'rebeccapurple',
        'Fibroblast': 'teal'
    }
elif specimen == 'human':
    colors = ['lightblue',"plum",'grey','mediumseagreen','gold','wheat','yellowgreen','olivedrab','rebeccapurple','teal','midnightblue','cyan','orangered','magenta']
    texts = ['EC',"PCT",'Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD']
    color_d = {
        'Podocyte' : 'midnightblue',
        'MC' : 'cyan',
        'EC' : 'lightblue',
        'PCT' : 'plum',
        'CD-IC' : 'olivedrab',
        'CD-PC' : 'yellowgreen',
        'Other_Immune' : 'grey',
        'Macrophage': 'mediumseagreen',
        'GC' : 'orangered',
        'MD' : 'magenta',
        'TAL' : 'gold',
        'DCT' : 'wheat',
        'vSMC': 'rebeccapurple',
        'Fibroblast': 'teal'
    }


ax.scatter(allcells_df['x'], allcells_df['y'], s=3,c=allcells_df['cell_type'].apply(lambda x: color_d[x]))
patches = [ plt.plot([],[], marker="o", ms=3, ls="", linewidth=0.6,color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
plt.legend(handles=patches, bbox_to_anchor=(1, 0.5), 
        loc='center left', ncol=1,fontsize='medium')
plt.setp(plt.gca().get_legend().get_texts(), fontsize='20',fontname='Arial')
plt.axis('off')
plt.savefig(out_path,dpi=300,bbox_inches='tight')
plt.close('all')