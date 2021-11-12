import pandas as pd
import numpy as np
import os
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
import argparse
np.random.seed(111)

parser = argparse.ArgumentParser()

parser.add_argument('--cell_type')
parser.add_argument('--section',help='cortex or medulla')

args = parser.parse_args()

cell_type = args.cell_type
section = args.section

# input_path is path to .txt file with list of unique ids of UMOD-KI arrays
UMODKI_array_ids = list(pd.read_csv(input_path,header=None)[0])
# input_path is path to .txt file with list of unique ids of UMOD-WT arrays
UMODWT_array_ids = list(pd.read_csv(input_path,header=None)[0])

array_ids = UMODKI_array_ids+UMODWT_array_ids

# for all UMOD-KI and UMOD-WT arrays,
# concatenate data matrices with sct-transformed UPR pathway gene counts for all curated beads of specified cell type in specified section

all_df = pd.DataFrame()
for array_id in array_ids:
    if array_id in UMODKI_array_ids:
        genotype = 'UMOD-KI'
    elif array_id in UMODWT_array_ids:
        genotype = 'UMOD-WT'
    
    # input_path is path to data matrix with beads x features for all curated cell types in array
    # features = {'barcode','x','y','cell_type','section'}
    input_path = '{}_allcells_df.csv'.format(array_id)
    allcells_df = pd.read_csv(input_path,index_col=0)
    
    df = allcells_df[allcells_df['section']==section].copy()
    df = df[df['cell_type']==cell_type].copy()
    
    input_path = '{array_id}_UPR_pathway_genes_sct_cts.feather'.format(array_id=array_id)
    sct_df = pd.read_feather(input_path)
    
    sct_df = sct_df.set_index(['genes'])
    sct_df = sct_df.T
    sct_df = sct_df.reset_index()
    sct_df = sct_df.rename(columns={'index':'barcode'})
    sct_df = df.merge(sct_df,on='barcode')
    sct_df['array_id'] = [array_id]*sct_df.shape[0]
    sct_df['genotype'] = [genotype]*sct_df.shape[0]
    sct_df = sct_df.set_index(['barcode','x','y','section','cell_type','array_id','genotype']).reset_index()
    
    all_df = pd.concat([all_df,sct_df])
all_df = all_df.reset_index()
all_df = all_df.drop(columns={'index'})

d = {'IRE1alpha-down': ['Yif1a','Tmem165'],
     'PERK-down': ['Hspe1','Rps26','Ppia','Ndufs5','Prdx1','Ptma','Dnaja1','Tubb4b','Uqcr11','Polr2l','Cox6b1','Rps10','Cox6a1'],
     'PERK-up': ['Mthfd2','Eif4ebp1','Xbp1','Ddit3','Trib3','Socs2','Cebpg','Eif1b','Eif1','Map1b','Gars','Pck2','Sesn2','Cth','Psph','Wars','Phgdh','Hax1','Psat1','Bex2','Lmo4','Fam89a','Tsc22d3','Rgs16','Map3k8','Idh1','Ccpg1','Pim1','Slc3a2','Snhg8'],
     'IRE1alpha-up': ['Dnajb9','Tmed2','Serp1','Vimp','Derl2','Slc35b1','Erlec1','Armcx3','Sec61a1','Sec61b','Ppib','Ssr2','Tmed9','Nans','Ostc','Ssr3','Ssr1'],
     'ATF6-up': ['Selk','Cdk2ap2','Hspa5','Herpud1','Sdf2l1','Dnajb11','Manf','Hsp90b1','Creld2','Pdia6','Pdia4','Calr','Dnajc3','Hyou1','Tmem50b'],
    }

order = ['Yif1a','Tmem165','Hspe1','Rps26','Ppia','Ndufs5','Prdx1','Ptma','Dnaja1','Tubb4b','Uqcr11','Polr2l','Cox6b1','Rps10','Cox6a1','Mthfd2','Eif4ebp1','Xbp1','Ddit3','Trib3','Socs2',
         'Cebpg','Eif1b','Eif1','Map1b','Gars','Pck2','Sesn2','Cth','Psph','Wars','Phgdh','Hax1','Psat1','Bex2','Lmo4','Fam89a','Tsc22d3','Rgs16','Map3k8','Idh1','Ccpg1','Pim1','Slc3a2',
         'Snhg8','Dnajb9','Tmed2','Serp1','Vimp','Derl2','Slc35b1','Erlec1','Armcx3','Sec61a1','Sec61b','Ppib','Ssr2','Tmed9','Nans','Ostc','Ssr3','Ssr1','Selk','Cdk2ap2','Hspa5',
         'Herpud1','Sdf2l1','Dnajb11','Manf','Hsp90b1','Creld2','Pdia6','Pdia4','Calr','Dnajc3',
         'Hyou1','Tmem50b']

array_to_mouse_d = {
    '191223_15':'1.1',
    '191223_17':'1.1',
    '191223_18':'1.1',
    '191223_19':'1.1',
    '191223_20':'1.1',
    '191223_21':'1.2',
    '191223_22':'1.2',
    '191223_23':'1.2',
    '191223_24':'1.2',
    '200102_09':'1.2',
    '200104_01':'2.1',
    '200104_02':'2.1',
    '200104_03':'2.1',
    '200104_04':'2.1',
    '200104_05':'2.1',
    '200104_06':'2.2',
    '200104_07':'2.2',
    '200104_09':'2.2',
    '200104_10':'2.2',
    '200104_14':'2.2',
    '200127_01':'3.1',
    '200127_05':'3.1',
    '200127_06':'3.2',
    '200127_07':'3.2',
    '200127_08':'3.2',
    '200127_09':'3.2',
    '200127_10':'3.2',
    '200131_13':'4.2',
    '200131_20':'4.2',
    '200131_21':'4.2',
    '200131_22':'4.2',
    '200131_23':'4.2',
    '200210_01':'5.2',
    '200210_02':'5.2',
    '200210_03':'5.2',
    '200210_04':'5.2'
}

temp=all_df.iloc[:,7:]
temp=temp.loc[:,(temp!=0).any(0)]
temp=temp.reset_index()
temp=temp.drop(columns={'index'})

info = all_df.iloc[:,0:7]
info=info.reset_index()
info=info.drop(columns={'index'})

dat = pd.DataFrame()
dat = pd.concat([info,temp],axis=1)

### track which arrays correspond to different mice and average gene expression by mouse
mouse_ids = []
for array_id in dat['array_id']:
    mouse_id = array_to_mouse_d[array_id]
    mouse_ids.append(mouse_id)
temp = dat.copy()
temp['mouse_id'] =  mouse_ids
temp=temp.drop(columns={'barcode','x','y','section','cell_type','array_id','genotype'})
mouse_avg = temp.groupby(['mouse_id']).mean()
genos = ['UMOD-WT','UMOD-KI','UMOD-WT','UMOD-KI','UMOD-WT','UMOD-KI','UMOD-KI','UMOD-KI']

mouse_avg['geno']=genos
mouse_avg=mouse_avg.set_index(['geno']).reset_index()

### DE analysis between UMOD-KI and UMOD-WT mice
df = mouse_avg.copy()
metadata = pd.DataFrame(df['geno'])
metadata['geno'] = metadata['geno'].astype('str')
metadata['geno'] = metadata['geno'].astype('category')

metadata=metadata.reset_index()
metadata=metadata.drop(columns={'index'})

counts_dat = df.iloc[:,1:]
counts_dat=counts_dat.reset_index()
counts_dat=counts_dat.drop(columns={'index'})

adata = sc.AnnData(X = counts_dat, obs = metadata)

n_genes = 74
sc.tl.rank_genes_groups(adata, groupby='geno', use_raw=True, 
                        method='wilcoxon', n_genes=n_genes)

cluster_labs = pd.DataFrame()
for i in ['UMOD-KI','UMOD-WT']:
    cluster_lab = pd.DataFrame([i]*n_genes)
    cluster_labs = pd.concat([cluster_labs,cluster_lab])

pvals = []
pvals_adj
gene_names = []
log_fcs = []
for i in range(n_genes):
    pvals.append(np.array([x for x in adata.uns['rank_genes_groups']['pvals'][i]]))
    pvals_adj.append(np.array([x for x in adata.uns['rank_genes_groups']['pvals_adj'][i]]))
    gene_names.append(np.array([x for x in adata.uns['rank_genes_groups']['names'][i]]))
    log_fcs.append(np.array([x for x in adata.uns['rank_genes_groups']['logfoldchanges'][i]]))
pval_rows = pd.DataFrame(np.array(pval_rows).T.ravel())
gene_names = pd.DataFrame(np.array(gene_names).T.ravel())
logs = pd.DataFrame(np.array(logs).T.ravel())

DE_dat = pd.DataFrame()
DE_dat['genes'] = np.array(gene_names[0])
DE_dat['pvals'] = np.array(pvals[0])
DE_dat['pvals_adj'] = np.array(pvals_adj[0])
DE_dat['cluster'] = np.array(cluster_labs[0])
DE_dat['log_fc'] = np.array(log_fcs[0])

# out_path is path to output file
out_path = '{cell_type}_{section}_UPR_pathway_DE_results.csv'.format(cell_type=cell_type,section=section)
DE_dat.to_csv(out_path)

DE_dat=DE_dat[DE_dat['pvals_adj']<0.05]
DE_dat=DE_dat[DE_dat['logs']>0]

### Order genes based on log foldchange, within each pathway
temp = [x for x in order if x not in list(DE_dat['genes'])]
nonsig_genes = pd.DataFrame()
nonsig_genes['genes'] = temp
nonsig_genes['pvals'] = [1.0]*nonsig_genes.shape[0]
nonsig_genes['pvals_adj'] = [1.0]*nonsig_genes.shape[0]
nonsig_genes['cluster'] = [np.nan]*nonsig_genes.shape[0]
nonsig_genes['log_fc'] = [0.0]*nonsig_genes.shape[0]

DE_dat = pd.concat([DE_dat,nonsig_genes])

pathways = []
for gene in DE_dat['genes']:
    for x,y in d.items():
        if gene in y:
            pathway = x
            break
    pathways.append(x)

DE_dat['pathway'] = pathways

DE_dat=DE_dat.sort_values(by=['log_fc'],ascending=False).groupby('pathway').head(DE_dat.shape[0])

### More reordering based on pathway and genotype
pathway_order = ['IRE1alpha-down','PERK-down', 'ATF6-up', 'IRE1alpha-up', 'PERK-up']
DE_dat_ord = pd.DataFrame()
for pathway in pathway_order:
    temp = DE_dat[DE_dat['pathway']==pathway].copy()
    temp1 = temp[temp['cluster']=='UMOD-KI']
    temp2 = temp[temp['cluster']=='UMOD-WT']
    temp3 = temp[temp['cluster'].isnull().values]
    temp = pd.concat([temp1,temp2,temp3])
    DE_dat_ord = pd.concat([DE_dat_ord,temp])

DE_dat_ord = DE_dat_ord.reset_index()
DE_dat_ord = DE_dat_ord.drop(columns={'index'})

DE_dat=DE_dat_ord[DE_dat_ord['pvals_adj']<0.05]
DE_dat=DE_dat_ord[DE_dat_ord['logs']>0]

sig_genes=list(DE_dat['genes'])
sig_gene_d = {}
for g in sig_genes:
    sig_gene_d[g] = '* '+g

order = ['geno']+list(DE_dat_ord['genes'])
mouse_avg=mouse_avg[order]
mouse_avg=mouse_avg.rename(columns={'geno':'genotype'})
temp1 = mouse_avg[mouse_avg['genotype']=='UMOD-WT']
temp2 = mouse_avg[mouse_avg['genotype']=='UMOD-KI']
mouse_avg = pd.concat([temp1,temp2])

genesonly=mouse_avg.iloc[:,1:]

pathways = []
for gene in genesonly.columns:
    for x,y in d.items():
        if gene in y:
            pathway = x
            break
    pathways.append(x)

pathways=pd.DataFrame(pathways)
pathways=pathways.rename(columns={0:'pathway'})
pathways=pathways.set_index(genesonly.columns)

pathways=pathways.rename(index=sig_gene_d)
genesonly=genesonly.rename(columns=sig_gene_d)

plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = 22
row_d = {'UMOD-KI':'firebrick','UMOD-WT':'dodgerblue'}
phenos = mouse_avg['genotype']
row_colors = phenos.map(row_d)
col_d = {'IRE1alpha-down':'orange','PERK-down': 'magenta','PERK-up':'forestgreen','IRE1alpha-up':'blueviolet','ATF6-up':'mediumseagreen'}
pathways = pathways['pathway']
col_colors = pathways.map(col_d)
cmap = sns.diverging_palette(260, 12,as_cmap=True)
sns.clustermap(genesonly,row_cluster=False,col_cluster=False,cmap='magma',xticklabels=1,yticklabels=False,row_colors=row_colors,col_colors=col_colors,cbar_pos=(0.05, 0.4, 0.05, 0.18),figsize=(27,20),standard_scale=1)
colors = ['dodgerblue',"firebrick"]
texts = ['WT','UMOD_KI']
patches = [ plt.plot([],[], marker="s", ms=10,ls='none', linewidth=0.6,color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
leg1=plt.legend(handles=patches, loc='best', bbox_to_anchor=(1.25, 1.9, 0.5, 0.5), framealpha=1,frameon=False)
colors = ["orange",'magenta','mediumseagreen','blueviolet','forestgreen']
texts = ['IRE1alpha-down','PERK-down','ATF6-up','IRE1alpha-up','PERK-up']
patches = [ plt.plot([],[], marker="s", ms=10,ls='none', linewidth=0.6,color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
leg2=plt.legend(handles=patches, loc='best', bbox_to_anchor=(1.85, 1.5, 0.5, 0.5), framealpha=1,frameon=False)
plt.gca().add_artist(leg1)

# out_path is path to output figure
out_path = '{cell_type}_{section}_avg_by_mouse_clustermap_magma.pdf'.format(cell_type=cell_type,section=section)
plt.savefig(out_path,dpi=300)
