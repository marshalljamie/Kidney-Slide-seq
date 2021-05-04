import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import path
import argparse
from sklearn.cluster import KMeans
from sklearn import metrics
import seaborn as sns

parser = argparse.ArgumentParser()

parser.add_argument('--array_id',help='unique id for array to be analyzed')
parser.add_argument('--cell_type')
parser.add_argument('--section',help='cortex or medulla')

args = parser.parse_args()

array_id = args.array_id
cell_type = args.cell_type
section = args.section

# input path is path to file with beads x features for the curated cell type of interest
# features = {'barcode','x','y'}
input_path = '{array_id}_{cell_type}_{section}_info.csv'.format(array_id=array_id,cell_type=cell_type,section=section)
if path.exists(input_path):
    celltype_info = pd.read_csv(input_path,index_col=0)

    X = np.array(celltype_info[['x','y']])

    # generate silhouette scores for kmeans, k in [2,100)
    sil_scores = []
    range_limit = celltype_info.shape[0]
    range_max = 50
    if (range_limit-1) > (range_max-1):
        range_limit = range_max
    print(range_limit)
    for i in range(2,range_limit):
        kmeans = KMeans(n_clusters=i, random_state=0).fit(X)
        cluster_centers = kmeans.cluster_centers_
        cluster_labels = kmeans.labels_
        sil_scores.append(metrics.silhouette_score(X, cluster_labels, metric='euclidean'))

    # find optimal k that maximizes silhouette score
    cluster_list=list(range(2,range_limit))
    n_clusters = cluster_list[sil_scores.index(max(sil_scores))]

    # kmeans with optimal k
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
    cluster_centers = kmeans.cluster_centers_
    cluster_labels = kmeans.labels_
    celltype_info['clusters'] = cluster_labels

    # out_path is path to output file
    out_path = '{array_id}_{cell_type}_{section}_info_clustered.csv'.format(array_id=array_id,cell_type=cell_type,section=section)
    celltype_info.to_csv(out_path)

    # visualize partitions
    d = {}
    for cluster in cluster_labels:
        d[cluster] = (cluster_centers[cluster][0],cluster_centers[cluster][1])

    # out_path is path to output figure
    out_path = '{array_id}_{cell_type}_{section}_partitions.png'.format(cell_type=cell_type,section=section,puck_id=puck_id)
    sns.set(rc={'figure.figsize':(5,5)})
    sns.set_style("white")
    p=sns.scatterplot(x=celltype_info.x, y=celltype_info.y, data=celltype_info, hue=celltype_info.clusters, palette='Spectral',alpha=0.3,ec='black',s=40,legend=False)
    sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
    p.set(xticks=[],yticks=[],xlabel=None,ylabel=None)
    plt.xlim(0,6000)
    plt.ylim(0,6000)
    plt.savefig(out_path,dpi=300,bbox_inches='tight')
    





