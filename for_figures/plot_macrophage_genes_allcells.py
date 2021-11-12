import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from os import path
np.random.seed(111)

parser = argparse.ArgumentParser()

parser.add_argument('--specimen',help='mouse or human')
parser.add_argument('--array_id',help='Unique id for array to be processed')

args = parser.parse_args()

specimen = args.specimen
array_id = args.array_id

# input_path is path to file with beads x features for all curated cell type calls in an array
# features = {'barcode','x','y','cell_type','section'}
input_path = '{array_id}_allcells.csv'.format(array_id=array_id)
allcells_df = pd.read_csv(input_path,index_col=0)
allcells_df = allcells_df.reset_index()
allcells_df = allcells_df.drop(columns={'index'})
allcells_df = allcells_df[['barcode','x','y','cell_type']].copy()
if specimen ==  'mouse':
    allcells_df = allcells_df.replace(['Bcell','NKT','DC'],'Other_Immune')
    allcells_df = allcells_df.replace({'Mesangial':'MC','Endothelial':'EC','CDPC':'CD-PC','CDIC':'CD-IC','Ren1':'GC','PCT1':'PCT_1','PCT2':'PCT_2'})
elif specimen == 'human':
    allcells_df = allcells_df.replace({'Mesangial':'MC','Endothelial':'EC','CDPC':'CD-PC','CDIC':'CD-IC','Ren1':'GC','Immune':'Other_Immune'})


macrophage_barcodes = list(allcells_df[allcells_df['cell_type'] == 'Macrophage']['barcode'])
#input_path is path to file with beads x features for all trem2-expressing macrophages in medulla
# features = {'barcode','x','y'}
features = 
input_path = '{array_id}_immune_Trem2_medulla_info.csv'.format(array_id=array_id)
macrophage_trem2_m = pd.DataFrame(columns=['barcode','x','y','cell_type'])
macrophage_trem2_c = pd.DataFrame(columns=['barcode','x','y','cell_type'])
if path.exists(input_path):
    macrophage_trem2_m = pd.read_csv(input_path,index_col=0)
# input_path is path to file with beads x features for all trem2-expressing macrophages in cortex
# features = {'barcode','x','y'}
input_path = '{array_id}_immune_Trem2_cortex_info.csv'.format(array_id=array_id)
if path.exists(input_path):
    macrophage_trem2_c = pd.read_csv(input_path,index_col=0)
macrophage_trem2 = pd.concat([macrophage_trem2_m,macrophage_trem2_c],sort=False)
macrophage_trem2['cell_type'] = ['Macrophage_Trem2+']*macrophage_trem2.shape[0]
macrophage_trem2 = macrophage_trem2[macrophage_trem2['barcode'].isin(macrophage_barcodes)].copy()
if specimen == 'human':
    # input_path is path to file with beads x features for all lyve1-expressing macrophages in cortex
    # features = {'barcode','x','y'}
    macrophage_lyve1_m = pd.DataFrame(columns=['barcode','x','y','cell_type'])
    macrophage_lyve1_c = pd.DataFrame(columns=['barcode','x','y','cell_type'])
    input_path = '{array_id}_immune_C1q_Lyve1_cortex_info.csv'.format(array_id=array_id)
    if path.exists(input_path):
        macrophage_lyve1_c = pd.read_csv(input_path)
    # input_path is path to file with beads x features for all lyve1-expressing macrophages in medulla
    # features = {'barcode','x','y'}
    input_path = '{array_id}_immune_C1q_Lyve1_medulla_info.csv'.format(array_id=array_id)
    if path.exists(input_path):
        macrophage_lyve1_m = pd.read_csv(input_path)
    macrophage_lyve1 = pd.concat([macrophage_trem2_m,macrophage_trem2_c],sort=False)
    macrophage_lyve1['cell_type'] = ['Macrophage_Lyve1+']*macrophage_lyve1.shape[0]
    macrophage_lyve1 = macrophage_lyve1[macrophage_lyve1['barcode'].isin(macrophage_barcodes)].copy()
    

if  specimen == 'mouse':
    allcells_df = pd.concat([allcells_df,macrophage_trem2],sort=False)
elif specimen == 'human':
    allcells_df = pd.concat([allcells_df,macrophage_trem2,macrophage_lyve1],sort=False)
allcells_df = allcells_df[['barcode','x','y','cell_type']].copy()
allcells_df = allcells_df.reset_index()
allcells_df = allcells_df.drop(columns={'index'})

if specimen == 'mouse':
    allcells_df['cell_type'] = pd.Categorical(allcells_df['cell_type'], ['EC','PCT_1','PCT_2','Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD','Macrophage_Trem2+'])
elif specimen == 'human':
    allcells_df['cell_type'] = pd.Categorical(allcells_df['cell_type'], ['EC','PCT','Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD','Macrophage_Trem2+','Macrophage_Lyve1+'])
allcells_df = allcells_df.sort_values("cell_type")

out_path = '{array_id}_allcells_macrophage_genes.pdf'.format(array_id=array_id)
subplot_kw = dict(xlim=(0, 6000), ylim=(0, 6000), autoscale_on=False)
fig, ax = plt.subplots(subplot_kw=subplot_kw,figsize=(10,10))

if specimen == 'mouse':
    colors = ['lightblue',"plum","pink",'grey','mediumseagreen','gold','wheat','yellowgreen','olivedrab','rebeccapurple','teal','midnightblue','cyan','orangered','magenta','red']
    texts = ['EC',"PCT_1","PCT_2",'Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD','Macrophage_Trem2+']
    sizes = [3]*(len(texts)-1)+[7]
    outline = colors[:-1]+['black']
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
        'Fibroblast': 'teal',
        'Macrophage_Trem2+': 'red'
    }
elif specimen == 'human':
    colors = ['lightblue',"plum",'grey','mediumseagreen','gold','wheat','yellowgreen','olivedrab','rebeccapurple','teal','midnightblue','cyan','orangered','magenta','red','red']
    texts = ['EC',"PCT",'Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD','Macrophage_Trem2+','Macrophage_Lyve1+']
    sizes = [3]*(len(texts)-2)+[7,7]
    outline = colors[:-2]+['black','black']
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
        'Fibroblast': 'teal',
        'Macrophage_Trem2+': 'red',
        'Macrophage_Lyve1+': 'red'
    }


ax.scatter(allcells_df['x'], allcells_df['y'], s=3,c=allcells_df['cell_type'].apply(lambda x: color_d[x]))
ax.scatter(macrophage_trem2['x'],macrophage_trem2['y'],s=50,c='red',ec='black',linewidth=0.6)
if specimen == 'human':
    ax.scatter(macrophage_lyve1['x'],macrophage_lyve1['y'],s=50,c='r',ec='black',linewidth=0.6)
patches = [ plt.plot([],[], marker="o", ms=sizes[i], ls="", mec=outline[i], linewidth=0.6,color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
plt.legend(handles=patches, bbox_to_anchor=(1, 0.5), 
        loc='center left', ncol=1)
plt.setp(plt.gca().get_legend().get_texts(), fontsize='20',fontname='Arial')
plt.axis('off')
plt.savefig(out_path,dpi=300,bbox_inches='tight')

if specimen == 'human':
    out_path = '{array_id}_allcells_macrophage_Lyve1.pdf'.format(array_id=array_id))
    subplot_kw = dict(xlim=(0, 6000), ylim=(0, 6000), autoscale_on=False)
    fig, ax = plt.subplots(subplot_kw=subplot_kw,figsize=(10,10))

    temp = allcells_df[~allcells_df['cell_type'].isin(['Macrophage_Trem2+'])]
    colors = ['lightblue',"plum",'grey','mediumseagreen','gold','wheat','yellowgreen','olivedrab','rebeccapurple','teal','midnightblue','cyan','orangered','magenta','red']
    texts = ['EC',"PCT",'Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD','Macrophage_Lyve1+']
    sizes = [3]*(len(texts)-1)+[7]
    outline = colors[:-1]+['black']
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
        'Fibroblast': 'teal',
        'Macrophage_Trem2+': 'red',
        'Macrophage_Lyve1+': 'red'
    }

    ax.scatter(temp['x'], temp['y'], s=3,c=temp['cell_type'].apply(lambda x: color_d[x]))
    ax.scatter(macrophage_lyve1['x'],macrophage_lyve1['y'],s=50,c='r',ec='black',linewidth=0.6)
    patches = [ plt.plot([],[], marker="o", ms=sizes[i], ls="", mec=outline[i], linewidth=0.6,color=colors[i], 
                    label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1, 0.5), 
            loc='center left', ncol=1)
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='20',fontname='Arial')
    plt.axis('off')
    plt.savefig(out_path,dpi=300,bbox_inches='tight')

    out_path = '{array_id}_allcells_macrophage_Trem2.pdf'.format(array_id=array_id)
    subplot_kw = dict(xlim=(0, 6000), ylim=(0, 6000), autoscale_on=False)
    fig, ax = plt.subplots(subplot_kw=subplot_kw,figsize=(10,10))

    temp = allcells_df[~allcells_df['cell_type'].isin(['Macrophage_Lyve1+'])]
    colors = ['lightblue',"plum",'grey','mediumseagreen','gold','wheat','yellowgreen','olivedrab','rebeccapurple','teal','midnightblue','cyan','orangered','magenta','red']
    texts = ['EC',"PCT",'Other_Immune','Macrophage','TAL','DCT','CD-PC','CD-IC','vSMC','Fibroblast','Podocyte','MC','GC','MD','Macrophage_Trem2+']
    sizes = [3]*(len(texts)-1)+[7]
    outline = colors[:-1]+['black']
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
        'Fibroblast': 'teal',
        'Macrophage_Trem2+': 'red',
        'Macrophage_Lyve1+': 'red'
    }

    ax.scatter(temp['x'], temp['y'], s=3,c=temp['cell_type'].apply(lambda x: color_d[x]))
    ax.scatter(macrophage_trem2['x'],macrophage_trem2['y'],s=50,c='r',ec='black',linewidth=0.6)
    patches = [ plt.plot([],[], marker="o", ms=sizes[i], ls="", mec=outline[i], linewidth=0.6,color=colors[i], 
                    label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1, 0.5), 
            loc='center left', ncol=1)
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='20',fontname='Arial')
    plt.axis('off')
    plt.savefig(out_path,dpi=300,bbox_inches='tight')
