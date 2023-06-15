import marsilea as ma
import marsilea.plotter as mp
import numpy as np
from anndata import AnnData
from scipy.sparse import issparse
from natsort import natsorted


class MarsileaHeatmap:

    def __init__(self,
                 adata: AnnData,
                 layer=None,
                 **kwds):
        if layer is None:
            cluster_data = adata.X
        else:
            cluster_data = adata.layers[layer]

        if issparse(cluster_data):
            cluster_data = cluster_data.toarray()

        self.adata = adata
        self.heatmap = ma.Heatmap(cluster_data, **kwds)

    PLOTTER_POOL = {
        "bar": mp.Numbers,
        "text": mp.Labels,
        "annot": mp.AnnoLabels,
        "colors": mp.Colors,
    }

    def _add_plot(self, slot, side, keys, plot, size=1, pad=.1, **kwds):
        data = getattr(self.adata, slot)[keys]
        plotter = self.PLOTTER_POOL[plot]
        self.heatmap.add_plot(side, plotter(data, **kwds), pad=pad, size=size)

    def add_left(self, keys, plot, size=1, pad=.1, **kwds):
        self._add_plot("obs", "left", keys, plot, size, pad, **kwds)

    def add_right(self, keys, plot, size=1, pad=.1, **kwds):
        self._add_plot("obs", "right", keys, plot, size, pad, **kwds)

    def add_top(self, keys, plot, size=1, pad=.1, **kwds):
        self._add_plot("var", "top", keys, plot, size, pad, **kwds)

    def add_bottom(self, keys, plot, size=1, pad=.1, **kwds):
        self._add_plot("var", "bottom", keys, plot, size, pad, **kwds)

    def _groupby(self, slot, keys, side, add_labels=True, **kwds):
        data = getattr(self.adata, slot)[keys]
        categorical_orders = natsorted(np.unique(data))
        self.heatmap.hsplit(labels=data, order=categorical_orders)
        if add_labels:
            self.heatmap.add_plot(side,
                                  mp.Chunk(categorical_orders, rotation=0, **kwds))

    def h_groupby(self, keys, add_labels=True, side="left", **kwds):
        self._groupby("obs", keys, side, add_labels, **kwds)

    def v_groupby(self, keys, add_labels=True, side="top", **kwds):
        self._groupby("var", keys, side, add_labels, **kwds)

    def add_dendrogram(self, side, **kwds):
        self.heatmap.add_dendrogram(side, **kwds)

    def add_title(self, *args, **kwds):
        self.heatmap.add_title(*args, **kwds)

    def legend(self, side="right", **kwds):
        self.heatmap.add_legends(side, **kwds)

    def show(self, return_axes=False):
        self.heatmap.render()
