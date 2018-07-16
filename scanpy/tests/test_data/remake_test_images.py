# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl

import scanpy.api as sc

adata = sc.datasets.krumsiek11()

# make heatmap
sc.pl.heatmap(adata, adata.var_names, 'cell_type', use_raw=False)
pl.savefig('master_heatmap.png', dpi=80, bbox_inches='tight')

# make heatmap with continues data
adata.obs['Gata2'] = adata.X[:, 0]
sc.pl.heatmap(adata, adata.var_names, 'Gata2', use_raw=False,
              num_categories=4, figsize=(4.5, 5))
pl.savefig('master_heatmap2.png', dpi=80, bbox_inches='tight')

# make dotplot
sc.pl.dotplot(adata, adata.var_names, 'cell_type', use_raw=False)
pl.savefig('master_dotplot.png', dpi=80, bbox_inches='tight')

# make dotplot with continuous data
sc.pl.dotplot(adata, adata.var_names, 'Gata2', use_raw=False,
              num_categories=7, figsize=(7, 2.5))
pl.savefig('master_dotplot2.png', dpi=80, bbox_inches='tight')

