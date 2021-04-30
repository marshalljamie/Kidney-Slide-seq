import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
import math

### Function defns

# get_markers_in_struct helper
# input: array of all coords and coords of thresholded marker
# returns one array with coordinates assigned marker and all other coordinates, with former assigned 1 and latter assigned 0 
def get_ingroup_outgroup(all_coords, marker_coords_thresholded):
    ingroup_coords = marker_coords_thresholded.tolist()
    all_coords = all_coords.tolist()
    outgroup_coords = [row for row in all_coords if row not in ingroup_coords]
    
    ingroup_coords = np.array(ingroup_coords, dtype='object')
    outgroup_coords = np.array(outgroup_coords, dtype='object')
    
    ingroup_coords_label = [(coord[0],coord[1],coord[2],1) for coord in ingroup_coords]
    outgroup_coords_label = [(coord[0],coord[1],coord[2],0) for coord in outgroup_coords]
    ingroup_coords_label.extend(outgroup_coords_label)
    coords_labeled = np.array(ingroup_coords_label,dtype='object')
    
    return(coords_labeled)

# erode helper
# input: knn adjacency matrix and array with coords labeled by group
# returns adjacency matrix where connections to outgroup coords along one axis have been removed
def remove_outgroup(nbrs, coords_labeled):
    outgroup_indices = np.where(coords_labeled[:,3]==0)[0]
    nbrs[outgroup_indices,:] = np.array([0.]*len(nbrs))
    
    return(nbrs)


# get_markers_in_struct helper
# input: knn adjacency matrix where extraneous edges have been removed, array with coords labeled by group,
# and threshold number of neighbors that coords in structure must have
# returns coordinates within structure 
def erode(nbrs_reduced, coords_labeled, threshold):
    # execute a version of erosion where only coordinates with a threshold number of neighbors are collected
    
    tot_nbrs = np.sum(nbrs_reduced,axis=1)
    indices = np.where(tot_nbrs > threshold)[0]
    
    eroded = coords_labeled[indices]
    indices_col = np.array([indices]).T
    eroded = np.append(eroded,indices_col,1)
    
    indices = indices.tolist()
    
    # retrieve edge points connected to coords collected in previous step 
    # (lost because they have connections to both outgroup and ingroup coords so they may not meet threshold)
    edge_pts = []
    for i in range(len(eroded)):
        i_nbrs = nbrs_reduced[int(eroded[i][4])] # nbrs or nbrs_reduced
        for j in range(len(i_nbrs)):
            if i_nbrs[j] == 1:
                if j not in indices:
                    index = np.array([j])
                    add = np.concatenate((coords_labeled[j],index))
                    indices.append(j)
                    edge_pts.append(add)
    
    edge_pts = np.array(edge_pts,dtype='object')
    eroded = np.concatenate((eroded, edge_pts))
    
    return(eroded)

# input: all coordinates, thresholded marker coordinates, desired number of neighbors for knn, 
# and threshold for number of neighbors a coord must have in order to be kept
# finds adjacency matrix with nearest neighbors of each coordinate, throws out unnecessary neighbors, and
# returns all coords of marker cell type contained within the desired structure
def get_markers_in_struct(all_coords, marker_coords_thresholded, n_neighbors, threshold):
    coords_labeled = get_ingroup_outgroup(all_coords, marker_coords_thresholded)
    coords_only = np.delete(coords_labeled, [0,3] ,1)
    
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree',n_jobs=-1).fit(coords_only)
    nbrs = nbrs.kneighbors_graph(coords_only).toarray()
    
    # remove unnecessary edges (loops and those connecting outgroup coords)
    nbrs_reduced = nbrs-np.identity(len(nbrs))
    nbrs_reduced = remove_outgroup(nbrs_reduced, coords_labeled)
    nbrs_reduced = remove_outgroup(nbrs_reduced.T, coords_labeled)
    
    eroded = erode(nbrs_reduced, coords_labeled, threshold)
    eroded = np.delete(eroded, [3,4], 1)

    return(eroded)    