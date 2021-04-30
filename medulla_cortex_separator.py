import numpy as np
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--input_dat')
parser.add_argument('--array_id',help='unique id for array to be analyzed')
parser.add_argument('--specimen',help='mouse or human')
parser.add_argument('--out_dir')

args = parser.parse_args()

# input_dat is path to file with beads x features for all curated cell type calls in an array
# features = {'barcode','x','y','cell_type'}
dat = args.dat
array_id = args.array_id
specimen = args.specimen
out_dir = args.out_dir

class SelectFromCollection:
    """
    Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax, collection, alpha_other=0.0):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    data = pd.read_csv(dat,index_col=0)
    data = np.array(data)
    barcodes = data[:,0]
    barcodes = [[x] for x in barcodes]
    barcodes = np.array(barcodes)
    labs = data[:,3]
    labs = [[x] for x in labs]
    labs = np.array(labs)

    subplot_kw = dict(xlim=(0, 6000), ylim=(0, 6000), autoscale_on=False)
    fig, ax = plt.subplots(subplot_kw=subplot_kw,figsize=(10,10))

    if specimen == 'mouse':
        colors = {
            'Podocyte' : 'midnightblue',
            'Mesangial' : 'cyan',
            'Endothelial' : 'lightblue',
            'PCT1' : 'plum',
            'PCT2' : 'pink',
            'CDIC' : 'olivedrab',
            'CDPC' : 'yellowgreen',
            'DC' : 'teal',
            'NKT': 'teal',
            'Bcell': 'teal',
            'Macrophage': 'mediumseagreen',
            'Ren1' : 'orangered',
            'MD' : 'magenta',
            'TAL' : 'gold',
            'DCT' : 'wheat',
            'vSMC': 'rebeccapurple',
            'Fibroblast': 'grey'
        }
    elif specimen == 'human':
        colors = {
            'Podocyte' : 'midnightblue',
            'Mesangial' : 'cyan',
            'Endothelial' : 'lightblue',
            'PCT' : 'plum',
            'CDIC' : 'olivedrab',
            'CDPC' : 'yellowgreen',
            'Immune' : 'teal',
            'Ren1' : 'orangered',
            'MD' : 'magenta',
            'TAL' : 'gold',
            'DCT' : 'wheat',
            'vSMC': 'rebeccapurple',
            'Fibroblast': 'grey'
        }

    color_lst = [colors[x[0]] for x in labs]

    pts = ax.scatter(data[:, 1], data[:, 2], s=3, c=color_lst)
    selector = SelectFromCollection(ax, pts)

    data = pd.DataFrame(data)
    data = data.rename(columns={0:'barcode',1:'x',2:'y',3:'cell_type'})

    def accept(event):
        if event.key == "enter":
        
            data_selected = np.array(selector.xys[selector.ind])
            labs_selected = labs[selector.ind]
            barcodes_selected = barcodes[selector.ind]
            data_selected = np.hstack((data_selected,labs_selected))
            data_selected = np.hstack((data_selected,barcodes_selected))
            data_selected = pd.DataFrame(data_selected)
            #data_selected = data_selected.rename(columns={0:'x',1:'y',2:'max_cell_type'})
            data_selected = data_selected.rename(columns={0:'x',1:'y',2:'cell_type',3:'barcode'})
            data_selected = data_selected[['barcode','x','y','cell_type']]
            out_path = os.path.join(out_dir,'{}_medulla_cells.csv'.format(array_id))
            data_selected.to_csv(out_path)

            out_path = os.path.join(out_dir,'{}_cortex_cells.csv'.format(array_id))
            data_not_selected = data[~data.index.isin(selector.ind)]
            data_not_selected.to_csv(out_path)
            
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()

    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Press enter to accept selected points.")

    plt.show()
