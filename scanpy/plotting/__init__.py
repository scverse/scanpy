from ._anndata import (
    scatter,
    violin,
    ranking,
    clustermap,
    tracksplot,
    dendrogram,
    correlation_matrix,
    heatmap,
)
from ._dotplot import DotPlot, dotplot
from ._matrixplot import MatrixPlot, matrixplot
from ._stacked_violin import StackedViolin, stacked_violin
from ._preprocessing import filter_genes_dispersion, highly_variable_genes

from ._tools.scatterplots import (
    embedding,
    pca,
    diffmap,
    draw_graph,
    tsne,
    umap,
    spatial,
)
from ._tools import pca_loadings, pca_scatter, pca_overview, pca_variance_ratio
from ._tools.paga import paga, paga_adjacency, paga_compare, paga_path
from ._tools import dpt_timeseries, dpt_groups_pseudotime
from ._tools import rank_genes_groups, rank_genes_groups_violin
from ._tools import (
    rank_genes_groups_dotplot,
    rank_genes_groups_heatmap,
    rank_genes_groups_stacked_violin,
    rank_genes_groups_matrixplot,
    rank_genes_groups_tracksplot,
)
from ._tools import sim
from ._tools import embedding_density

from ._rcmod import set_rcParams_scanpy, set_rcParams_defaults
from . import palettes

from ._utils import matrix
from ._utils import timeseries, timeseries_subplot, timeseries_as_heatmap

from ._qc import highest_expr_genes
