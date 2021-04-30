import pandas as pd
import numpy as np
import os
from IPython.display import display
import seaborn as sns
import scipy.stats
import scipy.optimize
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import NMF
from matplotlib import colors
import matplotlib.patches as mpatches
import collections
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse
import random
random.seed(111)
np.random.seed(111)

# Pipeline is from Rodriques et al., 2019
# Stripped of portions of code not relevant to cell type loadings
# Full source code found here: https://science.sciencemag.org/content/363/6434/1463
parser = argparse.ArgumentParser(description='NMF-reg pipeline for Slide-seq')

# matrix of format beads x features, where features = {'barcode','x','y'}
parser.add_argument('--ss_coords', help = 'Slide-seq array coordinates')
# matrix of format beads x features, where features = {'barcode','gene1','gene2',...'genen'}
parser.add_argument('--ss_cts', help = 'Slide-seq array umi matrix')
# matrix of format cells x features, where features = {'gene1','gene2',...,'genen'}
parser.add_argument('--sc_cts', help = 'Single-cell reference umi matrix')
# matrix of format cells x features, where features = {'clusters'}
parser.add_argument('--sc_clusters', help = 'Single-cell reference annotated clusters')
parser.add_argument('--specimen', help = "Specify \'human\' or \'mouse\'")
parser.add_argument('--k', help = 'Size of reduced-feature space', type=int)
parser.add_argument('--array_id', help = 'Unique id for array to-be-procesed')

args = parser.parse_args()
ss_coords = args.ss_coords
ss_cts = args.ss_cts
sc_cts = args.sc_cts
sc_clusters = args.sc_clusters
specimen = args.specimen
k = args.k
array_id = args.array_id

### Function definitions (complete set of functions from NMFreg developers; not all used for our purposes)
def plot_one_gene(gene):
    plt.figure(figsize=(10, 10))
    pyplot.set_cmap('viridis_r')
    plt.scatter(coords['x'], coords['y'], c=counts[gene], s=2, alpha=0.6)
    plt.axis('equal')
    plt.title('{}'.format(gene))
    plt.colorbar()
    plt.show()
    
def deconv_factor_to_celltype(row, adict, K, num_atlas_clusters):
    nc = num_atlas_clusters
    tmp_list = [0]*nc
    for key in range(K):
        item = adict[key] - 1
        tmp_list[item] += row[key]**2
    return pd.Series(np.sqrt(tmp_list))


def NMFreg(counts, coords, size, metacell_dict, gene_intersection, 
           num_atlas_clusters, celltype_to_factor_dict, 
           celltype_dict, plot_size_dict):

    puckcounts = counts[['barcode'] + gene_intersection]
    puckcounts = puckcounts.set_index(counts['barcode'])
    puckcounts = puckcounts.drop('barcode', axis=1)

    cell_totalUMI = np.sum(puckcounts, axis = 1)
    puckcounts_cellnorm = np.divide(puckcounts, cell_totalUMI[:,None])
    puckcounts_scaled = StandardScaler(with_mean=False).fit_transform(puckcounts_cellnorm)

    XsT = puckcounts_scaled.T

    Hs_hat = []
    for b in tqdm(range(XsT.shape[1])):
        h_hat = scipy.optimize.nnls(WaT, XsT[:, b])[0]
        if b == 0:
            Hs_hat = h_hat
        else:
            Hs_hat = np.vstack((Hs_hat, h_hat))

    Hs = pd.DataFrame(Hs_hat)
    Hs['barcode'] = puckcounts.index.tolist()

    Hs_norm = StandardScaler(with_mean=False).fit_transform(Hs.drop('barcode', 
                                                                    axis=1))

    Hs_norm = pd.DataFrame(Hs_norm)
    Hs_norm['barcode'] = puckcounts.index.tolist()

    
    maxloc_s = Hs_norm.drop('barcode', axis=1).values.argmax(axis=1)
    barcode_clusters = pd.DataFrame()
    barcode_clusters['barcode'] = Hs_norm['barcode']
    barcode_clusters['max_factor'] = maxloc_s

    barcode_clusters['atlas_cluster'] = barcode_clusters['barcode']

    for c in range(1, num_atlas_clusters + 1):
        condition = np.isin(barcode_clusters['max_factor'], 
                            celltype_to_factor_dict[c])
        barcode_clusters['atlas_cluster'][condition] = c       
        
    bead_deconv_df = Hs_norm.apply(lambda x: deconv_factor_to_celltype(row=x, 
                                            adict=factor_to_celltype_dict,
                                            K=K,
                                            num_atlas_clusters=num_atlas_clusters), 
                                   axis = 1)
    bead_deconv_df.insert(0, 'barcode', Hs_norm['barcode'])
    bead_deconv_df.columns = ['barcode'] + (bead_deconv_df.columns[1:]+1).tolist()
    bead_deconv_df = pd.DataFrame(bead_deconv_df)
    bead_deconv_df = bead_deconv_df.rename(columns = celltype_dict)
    
    maxloc_ct = bead_deconv_df.drop('barcode', axis=1).values.argmax(axis=1)+1
    bead_maxct_df = pd.DataFrame()
    bead_maxct_df['barcode'] = bead_deconv_df['barcode']
    bead_maxct_df['max_cell_type'] = maxloc_ct
    
    return Hs, Hs_norm, puckcounts, bead_deconv_df, barcode_clusters, bead_maxct_df, puckcounts_scaled

def deconv_factor_to_celltype_sum(row, adict, K, num_atlas_clusters):
    nc = num_atlas_clusters
    tmp_list = [0]*nc
    for key in range(K):
        item = adict[key] - 1
        tmp_list[item] += row[key]
    return pd.Series(tmp_list)

def deconv_factor_to_celltype_l2(row, adict, K, num_atlas_clusters):
    nc = num_atlas_clusters
    tmp_list = [0]*nc
    for key in range(K):
        item = adict[key] - 1
        tmp_list[item] += row[key]**2
    return pd.Series(np.sqrt(tmp_list))

def deconv_factor_to_celltype_mean(row, adict, K, num_atlas_clusters):
    nc = num_atlas_clusters
    tmp_list = [0]*nc
    for key in range(K):
        item = adict[key] - 1
        tmp_list[item] += row[key]
    num_fact = list(collections.OrderedDict(sorted(collections.Counter(adict.values()).items())).values()) 
    mean_tmp_list = np.divide(tmp_list, num_fact)
    return pd.Series(mean_tmp_list)

def cell_deconv(collapse):
    if(collapse=='l2'):
        tmp_df = Ha_norm.drop('cellname', axis=1).apply(lambda x: deconv_factor_to_celltype_l2(row=x, 
                                                            adict=factor_to_celltype_dict,
                                                            K=K,
                                                            num_atlas_clusters=num_atlas_clusters), 
                                                        axis = 1)
    
    if(collapse=='sum'):
        tmp_df = Ha_norm.drop('cellname', axis=1).apply(lambda x: deconv_factor_to_celltype_sum(row=x, 
                                                            adict=factor_to_celltype_dict,
                                                            K=K,
                                                            num_atlas_clusters=num_atlas_clusters), 
                                                        axis = 1)
    
    if(collapse=='mean'):
        tmp_df = Ha_norm.drop('cellname', axis=1).apply(lambda x: deconv_factor_to_celltype_mean(row=x, 
                                                            adict=factor_to_celltype_dict,
                                                            K=K,
                                                            num_atlas_clusters=num_atlas_clusters), 
                                                        axis = 1)
    
    tmp_df.insert(0, 'cellname', Ha_norm['cellname'])
    tmp_df.columns = ['cellname'] + (tmp_df.columns[1:]+1).tolist()
    tmp_df = pd.DataFrame(tmp_df)
    tmp_df = tmp_df.rename(columns = celltype_dict)

    maxloc_cellt = tmp_df.drop('cellname', axis=1).values.argmax(axis=1)+1
    cell_maxct_df = pd.DataFrame()
    cell_maxct_df['cellname'] = tmp_df['cellname']
    cell_maxct_df['max_cell_type'] = maxloc_cellt

    mismatch_df = cell_clusters[cell_maxct_df['max_cell_type'].values != 
                                cell_clusters['cluster'].values]
    #print('num mismatched: {}'.format(cell_clusters[cell_maxct_df['max_cell_type'].values != cell_clusters['cluster'].values].shape[0]))

    plt.figure(figsize=(4, 4))
    plt.hist(mismatch_df['cluster'])
    plt.close('all')

    return tmp_df, mismatch_df, cell_maxct_df

def plot_bar_cellt(cell_deconv_df_norm, cell_maxct_df, metacell_dict):
    for key, value in metacell_dict.items():
        ct_df = cell_deconv_df_norm[cell_maxct_df['max_cell_type']==int(key)]
        plt.figure(figsize=(4, 4))
        plt.bar(x=range(int(num_atlas_clusters)),
                height=np.sum(ct_df.drop(['cellname'], axis=1), axis=0),
                tick_label = list(metacell_dict.values()))
        plt.title(value)
        plt.xticks(rotation=90)
        plt.show()
        
def plot_hist_TF(cell_deconv_df_norm, cell_maxct_df, metacell_dict):
    posneg_dict = {}
    for key, value in metacell_dict.items():
        pos = cell_deconv_df_norm[value][cell_maxct_df['max_cell_type']==int(key)]
        neg = cell_deconv_df_norm[value][cell_maxct_df['max_cell_type']!=int(key)]
        posneg_dict[key] = [pos, neg]
        plt.figure(figsize=(4, 4))
        plt.hist(pos, range=(0, 1), color='green', alpha=0.6, density=True)
        plt.hist(neg, range=(0, 1), color='red', alpha=0.6, density=True)
        plt.title(value)
        plt.xticks(rotation=90)
        plt.close('all')
    return posneg_dict

def func_thresh_certainty(bead_deconv_df_norm, keep_thresh_df, metacell_dict):
    for key, value in metacell_dict.items():
        bool_df = keep_thresh_df[int(key)-1]
        ct_indx = list(bead_deconv_df_norm.index[bool_df.index])
        res = np.multiply(bead_deconv_df_norm['maxval'].loc[ct_indx], bool_df)
        bead_deconv_df_norm['thresh_ct'].loc[ct_indx] = res
        
    return bead_deconv_df_norm

def maxval_func(row): 
    return row[ct_names][metacell_dict[str(row['max_cell_type'])]]

def plot_boolean(size, coords, bead_maxct_df,
                 plot_size_dict, metacell_dict):
    for key, value in metacell_dict.items():
        boolcol = (bead_maxct_df['max_cell_type']==int(key))
        sub_df = bead_maxct_df.copy()
        sub_df['bool'] = boolcol

        plt.figure(figsize=(12, 12))
        plt.set_cmap('copper_r')
        plt.scatter(coords['x'], coords['y'], c=sub_df['bool'], 
                    s=plot_size_dict[size], alpha=0.6)
        plt.title('{} {}um'.format(value, size))
        plt.axis('equal')
        plt.show()
        
def plot_boolean_thresh(size, coords, bead_maxct_df, bead_deconv_df_norm,
                 plot_size_dict, metacell_dict):
    bool_col = bead_deconv_df_norm['thresh_ct']==0
    for key, value in metacell_dict.items():
        boolcol = bead_maxct_df['max_cell_type']==int(key)
        sub_df = bead_maxct_df.copy()
        sub_df['bool'] = boolcol
        bool_col_ct = np.multiply(bool_col, boolcol)

        plt.figure(figsize=(12, 12))
        plt.set_cmap('copper_r')
        plt.scatter(coords['x'], coords['y'], c=sub_df['bool'], 
                    s=plot_size_dict[size], alpha=0.6)
        plt.scatter(coords[bool_col_ct]['x'], 
                    coords[bool_col_ct]['y'], 
                    c='lightgray', s=plot_size_dict[size], alpha=1)
        plt.title('{} {}um'.format(value, size))
        plt.axis('equal')
        #plt.savefig(args.out+"/"+args.puck_id+"_"+'{}'.format(value)+"_thresh_binary.png",bbox_inches='tight',dpi=200)
        #plt.close('all')
        
def plot_ct_loadings(coords, size, puckcounts, bead_deconv_df, plot_size_dict):
    barcode_totalloading = np.sum(bead_deconv_df.drop('barcode', axis=1), 
                                  axis = 1)
    bead_deconv_df_norm = np.true_divide(bead_deconv_df.drop('barcode', axis=1), 
                                         barcode_totalloading[:,None])
    bead_deconv_df_norm['barcode'] = puckcounts.index.tolist()

    deconv_sub_df = bead_deconv_df_norm.drop('barcode', axis=1)

    for indx, col in deconv_sub_df.iteritems():
        plt.figure(figsize=(10, 10))
        plt.set_cmap('viridis_r')
        plt.scatter(coords['x'], coords['y'], c=bead_deconv_df_norm[indx], 
                    s=plot_size_dict[size], alpha=0.6)
        plt.title('{}'.format(indx))
        plt.axis('equal')
        plt.colorbar()
        plt.clim(0,1)
        #plt.show()
        #plt.savefig(args.out+"/"+args.puck_id+"_"+'{}'.format(indx)+"_loadings.png",bbox_inches='tight',dpi=200)
        #plt.close('all')
                                                                                            
    return bead_deconv_df_norm

def plot_allct(coords, size, bead_deconv_df_norm, 
               bead_maxct_df, plot_size_dict, num_atlas_clusters):
    bead_deconv_df_norm['max_cell_type'] = bead_maxct_df['max_cell_type']
    bead_deconv_df_norm['maxval'] = bead_deconv_df_norm.apply(maxval_func, 
                                                              axis=1)
    
    df_clust = pd.DataFrame(columns=['x','y','label'])
    df_clust['x'] = coords['x']
    df_clust['y'] = coords['y']
    df_clust['label'] = bead_deconv_df_norm['max_cell_type']

    facet = sns.lmplot(data=df_clust, x='x', y='y', hue='label', 
                       fit_reg=False, legend=False, legend_out=True,
                       palette = sns.color_palette("tab10", 
                                                   int(num_atlas_clusters)),
                       size = 10, scatter_kws={"s": 2*plot_size_dict[size]})
    #add a legend
    leg = facet.ax.legend(bbox_to_anchor=[1, 0.75],
                             title="label", fancybox=True)
    #change colors of labels
    for i, text in enumerate(leg.get_texts()):
        plt.setp(text, color = sns.color_palette("tab10", 
                                                 int(num_atlas_clusters))[i])
    return df_clust, bead_deconv_df_norm

def plot_certainty(coords, size, bead_deconv_df_norm, plot_size_dict):
    plt.figure(figsize=(12, 12))
    plt.set_cmap('Reds')
    plt.scatter(coords['x'], coords['y'], c=bead_deconv_df_norm['maxval'], 
                s=plot_size_dict[size], alpha=1)
    plt.title('Purity of most prevalent cell type per bead, {}um'.format(size))
    plt.axis('equal')
    plt.colorbar()
    plt.clim(0,1)
    plt.show()
    
def plot_certainty_thresh(coords, size, bead_deconv_df_norm, plot_size_dict):
    bool_col = bead_deconv_df_norm['thresh_ct']==0
    plt.figure(figsize=(12, 12))
    plt.set_cmap('Reds')
    plt.scatter(coords['x'], coords['y'], c=bead_deconv_df_norm['maxval'], 
                s=plot_size_dict[size], alpha=1)
    plt.colorbar()
    plt.scatter(coords[bool_col]['x'], 
                    coords[bool_col]['y'], 
                    c='lightgray', s=plot_size_dict[size], alpha=1)
    plt.title('Purity of most prevalent cell type per bead, {}um'.format(size))
    expr = round(100*np.divide(coords[bead_deconv_df_norm['thresh_ct']!=0].shape[0],
                               coords.shape[0]), 2)
    plt.xlabel('Single celltype beads: {}%'.format(expr))
    plt.axis('equal')
    plt.clim(0,1)
    plt.show()
    
def plot_certainty_perct(coords, size, bead_deconv_df_norm, df_clust, 
                         plot_size_dict, metacell_dict):
    for key, value in metacell_dict.items():
        bool_col = df_clust['label']==int(key)
        ct_df = bead_deconv_df_norm[bool_col]
        plt.figure(figsize=(12, 12))
        plt.set_cmap('Reds')
        plt.scatter(df_clust['x'], df_clust['y'], c='white', 
                    edgecolors='gray', linewidths=0.25, 
                    s=plot_size_dict[size], alpha=0.6)
        plt.scatter(coords[bool_col]['x'], coords[bool_col]['y'], 
                    c=ct_df['maxval'], s=plot_size_dict[size], alpha=0.9)
        plt.title('Purity of {} per bead, {}um'.format(value, size))
        plt.axis('equal')
        plt.colorbar()
        plt.clim(0,1)
        plt.show()
        
def plot_certainty_perct_thresh(coords, size, bead_deconv_df_norm, 
                                df_clust, bead_maxct_df, 
                                plot_size_dict, metacell_dict):
    keep_thresh_df = {}
    remove_thresh_df = {}
    bead_deconv_df_norm['max_cell_type'] = bead_maxct_df['max_cell_type']
    bead_deconv_df_norm['maxval'] = bead_deconv_df_norm.apply(maxval_func, 
                                                              axis=1)
    
    for key, value in metacell_dict.items():
        bool_col = (df_clust['label']==int(key)).tolist()
        ct_df = bead_deconv_df_norm[bool_col]
        ct_df['col'] = ct_df['maxval'].apply(lambda x: 0 if x <= thresh_certainty[int(key)-1] else x)
        keep_thresh_df[int(key)-1] = ct_df['maxval'] > thresh_certainty[int(key)-1]
        remove_thresh_df[int(key)-1] = ct_df['maxval'] <= thresh_certainty[int(key)-1]
        
        plt.figure(figsize=(12, 12))
        plt.set_cmap('Reds')
        plt.scatter(df_clust['x'], df_clust['y'], c='white', 
                    edgecolors='gray', linewidths=0.25, 
                    s=plot_size_dict[size], alpha=0.6)
        plt.scatter(coords[bool_col]['x'], coords[bool_col]['y'], 
                    c=ct_df['col'], s=plot_size_dict[size], alpha=0.9)
        plt.colorbar();
        plt.scatter(coords[bool_col]['x'][ct_df['col']==0], 
                    coords[bool_col]['y'][ct_df['col']==0], 
                    c='gray', s=plot_size_dict[size], alpha=0.4)
        
        plt.title('Purity of {} per bead, {}um'.format(value, size))
        plt.axis('equal')
        plt.clim(0,1)
        plt.close('all')
        #plt.show()
        
    return keep_thresh_df, remove_thresh_df

def plot_pie(size, bead_deconv_df_norm, keep_thresh_df, remove_thresh_df,
             thresh_certainty, metacell_dict):
    for key, value in metacell_dict.items():
        thresh_cert = thresh_certainty[int(key)-1]
        total = np.sum(bead_deconv_df_norm['max_cell_type']==int(key))
        keep = np.sum(keep_thresh_df[int(key)-1])
        remove = np.sum(remove_thresh_df[int(key)-1])
        print(value)
        print('total: {}'.format(total))
        print('passed certainty check: {}'.format(keep))
        print('did not pass: {}'.format(remove))
        print('certainty thresh: {}'.format(round(thresh_cert, 4)))
        print('____________')
        plt.figure(figsize=(6, 6))
        plt.pie(x=[keep, remove],
               labels = ['single celltype', 'mixed'],
               colors = ['xkcd:deep red', 'gray'],
               wedgeprops={'alpha':1},
               autopct =lambda x:'{:d}'.format(int(round(x*total/100))))
        title_str='{}um, {} purity threshold {}'.format(size, 
                                                        value, 
                                                        round(thresh_cert, 4))
        plt.title(title_str)
        plt.show()
        
def plot_bar_ct(bead_deconv_df_norm, df_clust, size, metacell_dict):
    for key, value in metacell_dict.items():
        ct_df = bead_deconv_df_norm[bead_deconv_df_norm['max_cell_type']==int(key)]
        plt.figure(figsize=(4, 4))
        plt.bar(x=range(int(num_atlas_clusters)),
                height=np.sum(ct_df.drop(['barcode', 'max_cell_type', 
                'maxval', 'thresh_ct'], axis=1), axis=0),
                tick_label = list(metacell_dict.values()))
        plt.title(value)
        plt.xticks(rotation=90)
        plt.show()
        
def plot_allct_thresh(coords, size, bead_deconv_df_norm, bead_maxct_df,
                      plot_size_dict, metacell_dict):
    bool_col = bead_deconv_df_norm['thresh_ct']==0
    plt.figure(figsize=(12, 12))
    for key, value in metacell_dict.items():
        boolcol = bead_maxct_df['max_cell_type']==int(key)
        #sub_df = bead_maxct_df.copy()
        #sub_df['bool'] = boolcol
        plt.scatter(coords[boolcol]['x'], coords[boolcol]['y'], 
                    #c=sns.color_palette("Paired", int(num_atlas_clusters))[int(key)-1], 
                    s=plot_size_dict[size], alpha=1)
        
        plt.scatter(coords[bool_col]['x'], 
                    coords[bool_col]['y'], 
                    c='lightgray', 
                    s=plot_size_dict[size], alpha=0.05)
    plt.title('All celltypes {}um'.format(size))
    plt.axis('equal')
    #plt.show()
    #plt.savefig(args.out+"/"+args.puck_id+"_allcells.png")
    #plt.close('all')
    return bead_deconv_df_norm

def get_df_clust(coords, bead_maxct_df):
    df_clust = pd.DataFrame(columns=['x','y','label'])
    df_clust['x'] = coords['x']
    df_clust['y'] = coords['y']
    df_clust['label'] = bead_maxct_df['max_cell_type']
    
    return df_clust

def plot_subset(user_genes,cell):
    color_list=['b','darkorange','green','r','darkviolet','saddlebrown','magenta','orchid','k','olive','teal','cornflowerblue','gold','lawngreen','lightsalmon','plum','peru','y','turquoise']
    c_num=0
    empty_beads=pcounts_and_coords[pcounts_and_coords[user_gene].sum(axis=1)==0]
    ax=empty_beads.plot(kind='scatter',x='xcoord',y='ycoord',color='lightgrey',alpha=.08,figsize=(10,10))
    for gene in user_genes:
        if gene in list(pcounts_and_coords.columns.values):
            tmp=pcounts_and_coords[pcounts_and_coords[gene]>0]
            if(len(tmp['xcoord'])!=0 and len(tmp['ycoord']) != 0):
                ax = tmp.plot(kind="scatter", x="xcoord",y="ycoord", color=color_list[c_num], label=gene,alpha=.5,ax=ax)
                c_num+=1
    plt.rcParams.update({'font.size': 20})
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-large')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(cell+" genes")
    #plt.savefig(args.out+"/"+args.puck_id+"_"+cell+"_genes.png",bbox_inches='tight',dpi=200)
    #plt.close('all')
    
### Preprocessing of file inputs
counts = pd.read_pickle(ss_cts)
coords = pd.read_csv(ss_coords)

bead_totals = np.sum(counts.drop(['barcode'], axis=1), axis=1)

atlas_dge = pd.read_csv(sc_cts)
cell_clusters = pd.read_csv(sc_clusters)

atlas_genes = atlas_dge.columns.tolist()

cell_clusters.columns = ['cellnames', 'cluster_name', 'cluster']

# mappings between cluster # and cell types
if (specimen).lower() == 'mouse':
    metacell_dict = {'1':"PCT-1",
                    '2':"Glom_endothelial",
                    '3':"Fenestrated_endothelial",
                    '4':"DCT",
                    '5':"Macrophages",
                    '6':"TAL",
                    '7':"PCT-2",
                    '8':"CD-B-IC",
                    '9':"Endothelial-1",
                    '10':"Podocyte",
                    '11':"CD-PC",
                    '12':"Mesangial",
                    '13':"NKT",
                    '14':"DC",
                    '15':"Renin_juxtaglomerular",
                    '16':"CD-A-IC",
                    '17':"Bcells_DC",
                    '18':"Endothelial-2",
                    '19':"PECs-1",
                    '20':"Thin_descending_limb"}

    celltype_dict = {1:"PCT-1",
                    2:"Glom_endothelial",
                    3:"Fenestrated_endothelial",
                    4:"DCT",
                    5:"Macrophages",
                    6:"TAL",
                    7:"PCT-2",
                    8:"CD-B-IC",
                    9:"Endothelial-1",
                    10:"Podocyte",
                    11:"CD-PC",
                    12:"Mesangial",
                    13:"NKT",
                    14:"DC",
                    15:"Renin_juxtaglomerular",
                    16:"CD-A-IC",
                    17:"Bcells_DC",
                    18:"Endothelial-2",
                    19:"PECs-1",
                    20:"Thin_descending_limb"}


# for human
elif (specimen).lower()=='human':
    metacell_dict = {'1':"CD-A-IC",
                    '2':"CD-B-IC",
                    '3':"CD-IC",
                    '4':"CD-PC",
                    '5':"CNT",
                    '6':"DCT",
                    '7':"Endothelial",
                    '8':"Fibroblasts-1",
                    '9':"IMCD",
                    '10':"Immune",
                    '11':"LOH",
                    '12':"Mito-rich",
                    '13':"PCT",
                    '14':"PEC",
                    '15':"Pericyte-2",
                    '16':"Pericyte-Mesangial",
                    '17':"Podocyte",
                    '18':"TAL",
                    '19':"tDL-1",
                    '20':"tDL-2",
                    '21':"Urothelial-1",
                    '22':"vSMC"}

    celltype_dict = {1:"CD-A-IC",
                    2:"CD-B-IC",
                    3:"CD-IC",
                    4:"CD-PC",
                    5:"CNT",
                    6:"DCT",
                    7:"Endothelial",
                    8:"Fibroblasts-1",
                    9:"IMCD",
                    10:"Immune",
                    11:"LOH",
                    12:"Mito-rich",
                    13:"PCT",
                    14:"PEC",
                    15:"Pericyte-2",
                    16:"Pericyte-Mesangial",
                    17:"Podocyte",
                    18:"TAL",
                    19:"tDL-1",
                    20:"tDL-2",
                    21:"Urothelial-1",
                    22:"vSMC"}

plot_size_dict = {10:4}
ct_names = list(metacell_dict.values())

gene_intersection = sorted(list(set(atlas_genes) & set(counts.columns.drop('barcode'))))
# remove genes that have zero counts in the sample
g_totals_sample = np.sum(counts.drop("barcode", axis=1)[gene_intersection], axis=0)
gene_intersection = list(g_totals_sample[g_totals_sample!=0].index)
atlasdge = atlas_dge[gene_intersection]

nonzero_beads = np.sum(counts[gene_intersection], axis=1)!=0
counts = counts[nonzero_beads]

K = args.k
random_state = 42

model = NMF(n_components=K, init='random', random_state = random_state)
Ha = model.fit_transform(atlasdge_scaled)
Wa = model.components_

Ha_norm = StandardScaler(with_mean=False).fit_transform(Ha)
Ha_norm = pd.DataFrame(Ha_norm)
Ha_norm['barcode'] = atlasdge.index.tolist()

maxloc = Ha_norm.drop('barcode', axis=1).values.argmax(axis=1)
cell_clusters['maxloc'] = maxloc

num_atlas_clusters = np.unique(cell_clusters['cluster']).size

factor_to_celltype_df = pd.DataFrame(0, index=range(1, num_atlas_clusters+1), 
                                     columns=range(K))
for k in range(K):
    n, bins, patches = plt.hist(cell_clusters['cluster'][cell_clusters['maxloc'] == k],
            range = (0.5, num_atlas_clusters+0.5), 
                                bins = int(num_atlas_clusters), 
                                facecolor='green', alpha=0.75)
    factor_to_celltype_df[k] = n.astype(int)


factor_to_celltype_df = factor_to_celltype_df.T

factor_total = np.sum(factor_to_celltype_df, axis = 1)
factor_to_celltype_df_norm = np.divide(factor_to_celltype_df, 
                                       factor_total[:,None])

maxloc_fc = factor_to_celltype_df.values.argmax(axis=1)
factor_to_celltype_dict = {factor : ctype + 1 for factor, ctype in enumerate(maxloc_fc)}

celltype_to_factor_dict = {}
for c in range(1, num_atlas_clusters + 1):
    celltype_to_factor_dict[c] = [k for k, v in factor_to_celltype_dict.items() if v == c]
    
WaT = Wa.T

Ha = pd.DataFrame(Ha)
Ha['cellname'] = atlasdge.index.tolist()

Ha_norm = StandardScaler(with_mean=False).fit_transform(Ha.drop('cellname', 
                                                                axis=1))

Ha_norm = pd.DataFrame(Ha_norm)
Ha_norm['cellname'] = atlasdge.index.tolist()

cell_deconv_df, mismatch_dfl2, cell_maxct_df = cell_deconv(collapse='l2')

cell_totalloading = np.sum(cell_deconv_df.drop('cellname', axis=1), axis = 1)
cell_deconv_df_norm = np.true_divide(cell_deconv_df.drop('cellname', axis=1), 
                                     cell_totalloading[:,None])
cell_deconv_df_norm['cellname'] = cell_deconv_df['cellname']

posneg_dict = plot_hist_TF(cell_deconv_df_norm=cell_deconv_df_norm, 
                           cell_maxct_df=cell_maxct_df, 
                           metacell_dict=metacell_dict)

thresh_certainty = [0]*num_atlas_clusters
for c in range(1, 1+num_atlas_clusters):
    thresh_certainty[c-1] = np.max(posneg_dict[str(c)][1])

Hs, Hs_norm, puckcounts, bead_deconv_df, barcode_clusters, bead_maxct_df, puckcounts_scaled = NMFreg(counts=counts, 
                                        coords=coords, size=10, 
                                        metacell_dict=metacell_dict, 
                                        gene_intersection=gene_intersection, 
                                        num_atlas_clusters=num_atlas_clusters, 
                                        celltype_to_factor_dict=celltype_to_factor_dict, 
                                        celltype_dict=celltype_dict, 
                                        plot_size_dict=plot_size_dict)

bead_deconv_df_norm = plot_ct_loadings(coords=coords, size=10, 
                        puckcounts=puckcounts, 
                        bead_deconv_df=bead_deconv_df,
                        plot_size_dict=plot_size_dict)

df_clust = get_df_clust(coords=coords, bead_maxct_df=bead_maxct_df)

keep_thresh_df, remove_thresh_df = plot_certainty_perct_thresh(coords=coords, 
                                      size=10, 
                                      bead_deconv_df_norm=bead_deconv_df_norm, 
                                      df_clust=df_clust,
                                      bead_maxct_df=bead_maxct_df, 
                                      plot_size_dict=plot_size_dict, 
                                      metacell_dict=metacell_dict)

df = pd.DataFrame(columns=['barcode', 'max_cell_type', 'x', 'y'])

for i in range(len(keep_thresh_df)):
    is_celltype_thresholded = keep_thresh_df[i]
    is_celltype = bead_maxct_df['max_cell_type'] == i + 1
    celltype_barcodes = bead_maxct_df[is_celltype]
    celltype_thresholded_barcodes = celltype_barcodes[is_celltype_thresholded]
    
    celltype_thresholded_coords = celltype_thresholded_barcodes.merge(coords, left_on='barcode', right_on='barcode')
    df = df.append(celltype_thresholded_coords)

bead_deconv_df_norm['thresh_ct'] = bead_deconv_df_norm['maxval']
bead_deconv_df_norm = func_thresh_certainty(bead_deconv_df_norm=bead_deconv_df_norm, 
                                            keep_thresh_df=keep_thresh_df,
                                            metacell_dict=metacell_dict)

bead_deconv_df_norm = plot_allct_thresh(coords=coords, size=10, 
                                    bead_deconv_df_norm=bead_deconv_df_norm, 
                                    bead_maxct_df=bead_maxct_df,
                                    plot_size_dict=plot_size_dict, 
                                    metacell_dict=metacell_dict)

#temp = bead_deconv_df_norm[['maxval','barcode']]
all_cell_loadings_thresholded = df.merge(bead_deconv_df_norm, on='barcode')
all_cell_loadings = coords.merge(bead_deconv_df_norm, on='barcode')

# out_path is path to output file
out_path = '{array_id}_thresholded_cell_loadings.csv'.format(array_id)
all_cell_loadings_thresholded.to_csv(out_path)
out_path = '{array_id}_all_cell_loadings.csv'.format(array_id)
all_cell_loadings.to_csv(out_path)
