import pandas as pd
import numpy as np
import os
import umap
import argparse
import hdbscan
import scanpy as sc
import matplotlib as mpl
import csv
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.sparse import csr_matrix
from sklearn.cluster import SpectralClustering
from os import path
np.random.seed(111)

parser = argparse.ArgumentParser()

parser.add_argument('--specimen',help='mouse or human')
parser.add_argument('--cell_type')
parser.add_argument('--section',help='cortex or medulla')
parser.add_argument('--data_form',help='bead or aggregate')
# only needs to be specified for mice
parser.add_argument('--genotype_list', nargs='+', help='list of genotypes to be compared, e.g. \'WT UMOD-KI\'')

# parameters for umap dimensionality reduction
# for descriptions on each of these parameters, see https://umap-learn.readthedocs.io/en/latest/parameters.html?highlight=parameters
parser.add_argument('--n_neighbors',type=int)
parser.add_argument('--min_dist',type=float)
parser.add_argument('--n_components',type=int)

# parameters for hdbscan clustering (clustering at aggregate level)
# for details, see https://umap-learn.readthedocs.io/en/latest/clustering.html
# can ignore if clustering at bead level
parser.add_argument('--min_samples',type=int) 
parser.add_argument('--min_cluster_size',type=int) 

args = parser.parse_args()

### inputs/inits
specimen = args.specimen
cell_type = args.cell_type
data_form = args.data_form
section = args.section
genotypes = args.genotype_list

# input path is path to data matrix with beads x features for curated beads (or aggregates) of specified cell type in specified section 
# belonging to arrays of diseased and healthy groups to-be-compared (i.e., (WT vs. UMOD_KI, BTBR wt/wt vs. BTBR ob/ob, etc.))
# for clustering on bead basis, features = {'barcode','array_id','genotype','batch','gene1','gene2',...,'genen'}; gene counts are sctransform residuals from Seurat
# for clustering on aggregate basis, features = {'array_id','genotype','batch','gene1','genen2',...,'genen'}; gene counts are averages of sctransform residuals from Seurat per structure
input_path = 'combined_cts.pkl'
combined_cts = pd.read_pickle(input_path)

### preprocessing
# set index to first column index that is a gene
index=0
combined_cts_genesOnly = combined_cts.iloc[:,index:]
combined_cts_arr = np.array(combined_cts_genesOnly)

### clustering
# set parameters for umap
n_neighbors = args.n_neighbors
min_dist = args.min_dist
n_components = args.n_components
# parameters for hdbscan
min_samples = args.min_samples
min_cluster_size = args.min_cluster_size
random_state = 42

# not always precisely reproducible when multiple cores requested (used for our purposes because vast amount of data on bead level across all arrays) 
# multiple runs will produce qualitatively reproducible results
standard_embedding = umap.UMAP(random_state=random_state).fit_transform(combined_cts_arr)
clusterable_embedding = umap.UMAP(
    n_neighbors=n_neighbors,
    min_dist=min_dist,
    n_components=n_components,
    random_state=random_state
).fit_transform(combined_cts_arr)

### spectral clustering if clustering on bead level
# set parameters for spectral clustering

if data_form == 'bead':
    n_clusters = 2
    assign_labels = 'discretize'
    random_state = 42

    labels = SpectralClustering(
        n_clusters=n_clusters,
        assign_labels=assign_labels,
        random_state=random_state).fit(clusterable_embedding)
    labels=labels.labels_
elif data_form == 'aggregate':
    labels = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size,
    ).fit_predict(clusterable_embedding)

clustered = (labels >= 0)
nonzero_clusters = np.unique(labels[clustered])
n_clusters = len(nonzero_clusters)

# output labels to out_path, path to file output
labels_out = pd.DataFrame(labels)
out_path = 'labels.csv'
labels_out.to_csv(out_path)

# output clusterable embedding, umap-transformed data to-be-clustered, to out_path, path to file output
clusterable_embedding_out = pd.DataFrame(clusterable_embedding)
clusterable_embedding_out = clusterable_embedding_out.rename(columns={0:'x',1:'y'})
out_path = 'clusterable_embedding.csv'
clusterable_embedding_out.to_csv(out_path)

# output standard embedding, data reduced to 2 dimensions for visualization purposes, to out_path, path to file output
out_path = 'standard_embedding.csv'
standard_embedding_out = pd.DataFrame(standard_embedding)
standard_embedding_out = standard_embedding_out.rename(columns={0:'x',1:'y'})
standard_embedding_out.to_csv(out_path)

# for plotting cluster labels in 2D umap space
plt.figure(figsize=(5,5))
colors = ['#1f77b4', '#ff7f0e','mediumseagreen','r','blueviolet','brown']
plt.scatter(standard_embedding[~clustered, 0],
            standard_embedding[~clustered, 1],
            c=(0.5,0.5,0.5),
            s=15,
            alpha=0.5)
for i in nonzero_clusters:
    plt.scatter(standard_embedding[(labels==i), 0],
                standard_embedding[(labels==i), 1],
                c=colors[i],
                ec = 'black',
                linewidth = 0.5,
                s=15)

plt.xlabel('UMAP-1',fontsize=12)
plt.ylabel('UMAP-2',fontsize=12)
plt.rcParams["axes.grid"] = False

colors = ['#1f77b4','#ff7f0e','lightgrey']
texts = ['cluster 0','cluster 1', 'not clustered']
outlines = ['black','black','grey']
patches = [ plt.plot([],[], marker="o", ms=10, mec=outlines[i],ls='none', linewidth=0.6,color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
plt.legend(handles=patches, loc='center left', framealpha=1,frameon=False,bbox_to_anchor=(1, 0.5),title='Cluster id')

# out_path is path to file output
out_path = 'clusters_umap.png'
plt.savefig(out_path,bbox_inches='tight',dpi=300)

### for plotting array disease state levels in 2D umap space
# for mice
if specimen == 'mouse':
    geno_cts = {}
    for geno in genotypes:
        ct = sum(combined_cts['genotype'][clustered]==geno.replace('-',''))
        geno_cts[geno] = ct
    print(geno_cts)

    print(geno_cts[genotypes[0]])
    plt.figure(figsize=(5,5))
    plt.scatter(standard_embedding[~clustered, 0],
                standard_embedding[~clustered, 1],
                c=(0.5,0.5,0.5),
                s=15,
                alpha=0.5)
    plt.scatter(standard_embedding[clustered,0][0:geno_cts[genotypes[0]]],
                standard_embedding[clustered,1][0:geno_cts[genotypes[0]]],
                c='r',
                ec='black',
                linewidth=0.5,
                s=15);
    plt.scatter(standard_embedding[clustered,0][geno_cts[genotypes[0]]:],
                standard_embedding[clustered,1][geno_cts[genotypes[0]]:],
                c='b',
                ec='black',
                linewidth=0.5,
                s=15);
    
    colors = ['b','r','lightgrey']
    texts = [genotypes[1],genotypes[0],'not clustered']
    outlines = ['black','black','grey']
    patches = [ plt.plot([],[], marker="o", ms=10, mec=outlines[i],ls='none', linewidth=0.6,color=colors[i], 
                    label="{:s}".format(texts[i]) )[0] for i in range(len(texts)) ]
    plt.legend(handles=patches, loc='center left', framealpha=1,frameon=False,bbox_to_anchor=(1, 0.5),title='Genotype')
    # out_path is path to file output
    out_path = 'disease_state_umap.png'
    plt.savefig(out_path,dpi=300,bbox_inches='tight')

    # compute proportion of each genotype in each cluster
    temp = combined_cts[['genotype']].copy()
    temp['labels'] = labels
    temp=temp[temp['labels']!=-1].copy()
    cluster_ids = np.unique(temp['clusteru'])
    genos = np.unique(temp['genotype'])
    result = []
    for cluster in cluster_ids:
        one_clust = temp[temp['labels']==cluster].copy()
        clust_tot = one_clust.shape[0]
        for geno in genos:
            one_geno=one_clust[one_clust['genotype']==geno].copy()
            geno_tot = one_geno.shape[0]
            prop = geno_tot/clust_tot
            result.append([geno,cluster,prop])
    # out_path is path to output file
    out_path = 'genotype_proportions_per_cluster.csv'
    result.to_csv(out_path)

elif specimen == 'human':
    # in our analysis, we sought to see where data points from DKD, ischemic, and normal humans distributed across clusters
    # mapping between unique array ids and disease status for humans
    arrayid_to_disease = {
        '200104_15':'unknown',
        '200104_16': 'unknown',
        '200104_17': 'unknown',
        '200104_18': 'unknown',
        '200104_19': 'DKD',
        '200104_20': 'DKD',
        '200104_21': 'DKD',
        '200104_23': 'DKD',
        '200115_01': 'unknown',
        '200115_02': 'unknown',
        '200115_03': 'unknown',
        '200115_04': 'unknown',
        '200115_05': 'unknown',
        '200115_06': 'unknown',
        '200115_07': 'unknown',
        '200115_09': 'unknown',
        '200115_10': 'unknown',
        '200115_11': 'unknown',
        '200115_12': 'unknown',
        '200115_14': 'unknown',
        '200115_15': 'Ischemic',
        '200115_16': 'Ischemic',
        '200115_17': 'Ischemic',
        '200115_18': 'Ischemic',
        '200113_07': 'unknown',
        '200113_08': 'unknown',
        '200113_09': 'unknown',
        '200113_10': 'unknown',
        '200113_11': 'Normal',
        '200113_12': 'Normal',
        '200121_01': 'Normal',
        '200121_03': 'Normal',
        '200131_24': 'unknown',
        '200131_25': 'unknown',
        '200131_26': 'unknown',
        '200205_13': 'unknown'
    }

    #puckids = list(combined_cts['array_id'])
    puckids = list(combined_cts['puck_id'])

    standard_embedding = pd.DataFrame(standard_embedding)
    standard_embedding = standard_embedding.rename(columns={0:'x',1:'y'})
    standard_embedding['disease_state'] = [arrayid_to_disease[x] for x in puckids]
    standard_embedding['disease_state'] = standard_embedding.disease_state.astype('category')
    standard_embedding['labels'] = labels

    unclustered = standard_embedding[~clustered]
    standard_embedding = standard_embedding[clustered]
    standard_embedding.to_csv('standard_embedding_mod.csv')

    plt.figure(figsize=(5,5))
    colors = {'Normal':'b','unknown':'black','DKD':'r','Ischemic':'purple'}
    plt.scatter(unclustered.x,unclustered.y,c='grey',alpha=0.5,s=15,ec='black',linewidth=0.5)
    plt.scatter(standard_embedding.x,standard_embedding.y,c=standard_embedding['disease_state'].map(colors),s=15,ec='black',linewidth=0.5)
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.xlabel('UMAP-1')
    plt.ylabel('UMAP-2')
    plt.rcParams["axes.grid"] = False

    texts = ['DKD','Ischemic','Normal','unknown']
    n = len(texts)
    texts = texts+['not clustered']
    outlines = ['black']*n+['grey']
    colors=['red','purple','b','black','grey']
    
    patches = [ plt.plot([],[], marker="o", ms=10, mec=outlines[i],ls='none', linewidth=0.6,color=colors[i], 
                    label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
    plt.legend(handles=patches, loc='center left', framealpha=1,frameon=False,bbox_to_anchor=(1, 0.5),title='Disease state')
    
    # out_path is path to output file
    out_path = 'disease_state_umap.png'
    plt.savefig(out_path,dpi=300,bbox_inches='tight')

    # compute proportion of each disease state in each cluster
    cluster_ids = np.unique(standard_embedding['labels'])
    disease_states = ['Normal','Ischemic','DKD','unknown']
    result = []
    for cluster in clusterids:
        one_clust = standard_embedding[standard_embedding['labels']==cluster].copy()
        clust_tot = one_clust.shape[0]
        for disease_state in disease_states:
            one_geno = one_clust[one_clust['disease_state']==disease_state].copy()
            geno_tot = one_geno.shape[0]
            prop = geno_tot/clust_tot
            result.append([disease_state,cluster,prop])
    # out_path is path to output file
    out_path = 'disease_state_proportions_per_cluster.csv'
    result.to_csv(out_path)

### Cluster DE
combined_cts['clusteru'] = labels
combined_cts = combined_cts.set_index(['clusteru']).reset_index()
counts_dat = combined_cts[combined_cts['clusteru']!=-1]
index = index+1
counts_dat = counts_dat.iloc[:,index:]

metadata = pd.DataFrame(combined_cts['clusteru'])
metadata = metadata[metadata['clusteru'] != -1]
metadata['clusteru'] = metadata['clusteru'].astype('str')
metadata['clusteru'] = metadata['clusteru'].astype('category')

metadata=metadata.reset_index()
metadata=metadata.drop(columns={'index'})
counts_dat=counts_dat.reset_index()
counts_dat=counts_dat.drop(columns={'index'})

adata = sc.AnnData(X = counts_dat, obs = metadata)

# out_path is path to output directory
out_path = ''
sc.settings.figdir = out_path

n_genes = 100
sc.tl.rank_genes_groups(adata, groupby='clusteru', use_raw=True, 
                        method='wilcoxon', n_genes=n_genes) # compute differential expression
sc.set_figure_params(fontsize=12,dpi_save=300,format='png')

# parses AnnData object for DE results (genes and p-values) 
cluster_labs = pd.DataFrame()
for i in nonzero_clusters:
    cluster_lab = pd.DataFrame([i]*n_genes)
    cluster_labs = pd.concat([cluster_labs,cluster_lab])

pvals = []
pvals_adj = []
gene_names = []
for i in range(n_genes):
    pvals.append(np.array([x for x in adata.uns['rank_genes_groups']['pvals'][i]]))
    pvals_adj.append(np.array([x for x in adata.uns['rank_genes_groups']['pvals_adj'][i]]))
    gene_names.append(np.array([x for x in adata.uns['rank_genes_groups']['names'][i]]))
pvals = pd.DataFrame(np.array(pvals).T.ravel())
pvals_adj = pd.DataFrame(np.array(pvals_adj).T.ravel())
gene_names = pd.DataFrame(np.array(gene_names).T.ravel())

DE_dat = pd.DataFrame()
DE_dat['genes'] = np.array(gene_names[0])
DE_dat['pvals'] = np.array(pvals[0])
DE_dat['pvals_adj'] = np.array(pvals_adj[0])
DE_dat['cluster'] = np.array(cluster_labs[0])

# out_path is path to output file
out_path = 'DE_dat.csv'
DE_dat.to_csv(out_path)

### save parameters used for dimensionality reduction and clustering
# for umap
# out_path is path to output file
out_path = 'umap_params.csv'
a_file = open(out_path, "w")
param_dict = {
            'n_neighbors' : n_neighbors,
            'min_dist' : min_dist,
            'n_components' : n_components,
            'random_state' : random_state,
        }
writer = csv.writer(a_file)

for key, value in param_dict.items():
    writer.writerow([key, value])
a_file.close()

if data_form == 'agg':
    # for hdbscan
    # out_path is path to output file
    out_path = 'hdbscan_params.csv'
    a_file = open(out_path, "w")
    param_dict = {
                'min_samples' : min_samples,
                'min_cluster_size' : min_cluster_size
            }
    writer = csv.writer(a_file)

    for key, value in param_dict.items():
        writer.writerow([key, value])
    a_file.close()
elif data_form == 'bead':
    # for spectral clustering
    # out_path is path to output file
    out_path = 'spectral_params.csv'
    a_file = open(out_path, "w")
    param_dict = {
            'n_clusters': n_clusters,
            'assign_labels': assign_labels,
            'random_state': random_state
    }
    writer = csv.writer(a_file)

    for key, value in param_dict.items():
        writer.writerow([key, value])
    a_file.close()
















