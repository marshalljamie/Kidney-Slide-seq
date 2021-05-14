import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--specimen',help='human or mouse')
parser.add_argument('--section',help='cortex or medulla')
parser.add_argument('--cell_type')

args = parser.parse_args()

section = args.section
cell_type = args.cell_type
specimen = args.specimen

if specimen == 'mouse':
    # .txt file with list of unique ids of DKD mouse arrays
    DKD_arrayids = list(pd.read_csv('mouse_DKD_arrayids.txt',header=None)[0])
    # .txt file with list of unique ids of WT mouse arrays
    WT_arrayids = list(pd.read_csv('mouse_WT_array_ids.txt',header=None)[0])
    # .txt file with list of unique ids of UMOD-KI arrays
    UMODKI_arrayids = list(pd.read_csv('mouse_UMODKI_arrayids.txt',header=None)[0])
    # .txt file with list of unique ids of UMOD-WT arrays
    UMODWT_arrayids = list(pd.read_csv('mouse_UMODWT_arrayids.txt',header=None)[0])
    array_ids = DKD_arrayids+WT_arrayids+UMODKI_arrayids+UMODWT_arrayids
elif specimen == 'human':
    # .txt file with list of unique ids of DKD human arrays
    DKD_arrayids = list(pd.read_csv('human_DKD_arrayids.txt',header=None)[0])
    # .txt file with list of unique ids of ischemic human arrays
    ischemic_arrayids = list(pd.read_csv('human_ischemic_arrayids.txt',header=None)[0])
    # .txt file with list of unique ids of healthy human arrays
    normal_arrayids = list(pd.read_csv('human_normal_arrayids.txt',header=None)[0])
    array_ids = DKD_arrayids+ischemic_arrayids+normal_arrayids


df_all = pd.DataFrame()
for array_id in array_ids:
    if specimen == 'mouse':
        if array_id in DKD_arrayids:
            genotype = 'DKD'
        elif array_id in WT_arrayids:
            genotype = 'WT'
        elif array_id in UMODKI_arrayids:
            genotype = 'UMOD-KI'
        elif array_id in UMODWT_arrayids:
            genotype = 'UMOD-WT'
    elif specimen == 'human':
        if array_id in DKD_arrayids:
            genotype = 'DKD'
        elif array_id in ischemic_arrayids:
            genotype = 'ischemic'
        elif array_id in normal_arrayids:
            genotype = 'normal'
    
    # input_path is path to file with beads x features for all curated cell type calls in an array
    # features = {'barcode','x','y','cell_type','section'}
    input_path = '{array_id}_allcells.csv'.format(array_id=array_id)
    allcells_info = pd.read_csv(input_path,index_col=0)
    allcells_info = allcells_info[allcells_info['section']==section].copy()

    celltype_info = celltype_info[celltype_info['cell_type']==cell_type].copy()
    num_celltype = celltype_info.shape[0]

    # input_path is path to file with all beads x features, with spatial outliers removed
    # features = {'barcode','x','y','section'}
    input_path = '{array_id}_coords_outliers_removed.csv'.format(array_id=array_id)
    coords_outliers_removed = pd.read_csv(input_path,index_col=0)
    coords_outliers_removed = coords_outliers_removed[coords_outliers_removed['section']==section].copy()
    tot = coords_outliers_removed.shape[0]
    
    pct = (num_celltype/tot)*100
    d = {'genotype':[genotype],'pct_{section}_{cell_type}'.format(section=section,cell_type=cell_type):[pct]}
    df = pd.DataFrame(d)
    df_all = pd.concat([df_all,df])

# out_path is path to output file
out_path = '{specimen}_{cell_type}_{section}_pct.csv'.format(specimen=specimen,cell_type=cell_type,section=section)
df_all.to_csv(out_path)