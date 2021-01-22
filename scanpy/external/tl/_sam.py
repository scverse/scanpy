"""\
Run the Self-Assembling Manifold algorithm
"""
from typing import Optional, Union, Tuple, Any

from anndata import AnnData

from ... import logging as logg
from ..._compat import Literal

try:
    from samalg import SAM
except ImportError:
    SAM = Any


def sam(
    adata: AnnData,
    max_iter: int = 10,
    num_norm_avg: int = 50,
    k: int = 20,
    distance: str = 'correlation',
    standardization: Literal['Normalizer', 'StandardScaler', 'None'] = 'StandardScaler',
    weight_pcs: bool = False,
    sparse_pca: bool = False,
    n_pcs: Optional[int] = 150,
    n_genes: Optional[int] = 3000,
    projection: Literal['umap', 'tsne', 'None'] = 'umap',
    inplace: bool = True,
    verbose: bool = True,
) -> Union[SAM, Tuple[SAM, AnnData]]:
    """\
    Self-Assembling Manifolds single-cell RNA sequencing analysis tool [Tarashansky19]_.

    SAM iteratively rescales the input gene expression matrix to emphasize
    genes that are spatially variable along the intrinsic manifold of the data.
    It outputs the gene weights, nearest neighbor matrix, and a 2D projection.

    The AnnData input should contain unstandardized, non-negative values.
    Preferably, the data should be log-normalized and no genes should be filtered out.


    Parameters
    ----------

    k
        The number of nearest neighbors to identify for each cell.

    distance
        The distance metric to use when identifying nearest neighbors.
        Can be any of the distance metrics supported by
        :func:`~scipy.spatial.distance.pdist`.

    max_iter
        The maximum number of iterations SAM will run.

    projection
        If 'tsne', generates a t-SNE embedding. If 'umap', generates a UMAP
        embedding. If 'None', no embedding will be generated.

    standardization
        If 'Normalizer', use sklearn.preprocessing.Normalizer, which
        normalizes expression data prior to PCA such that each cell has
        unit L2 norm. If 'StandardScaler', use
        sklearn.preprocessing.StandardScaler, which normalizes expression
        data prior to PCA such that each gene has zero mean and unit
        variance. Otherwise, do not normalize the expression data. We
        recommend using 'StandardScaler' for large datasets with many
        expected cell types and 'Normalizer' otherwise. If 'None', no
        transformation is applied.

    num_norm_avg
        The top 'num_norm_avg' dispersions are averaged to determine the
        normalization factor when calculating the weights. This prevents
        genes with large spatial dispersions from skewing the distribution
        of weights.

    weight_pcs
        If True, scale the principal components by their eigenvalues. In
        datasets with many expected cell types, setting this to False might
        improve the resolution as these cell types might be encoded by lower-
        variance principal components.

    sparse_pca
        If True, uses an implementation of PCA that accepts sparse inputs.
        This way, we no longer need a temporary dense copy of the sparse data.
        However, this implementation is slower and so is only worth using when
        memory constraints become noticeable.

    n_pcs
        Determines the number of top principal components selected at each
        iteration of the SAM algorithm. If None, this number is chosen
        automatically based on the size of the dataset. If weight_pcs is
        set to True, this parameter primarily affects the runtime of the SAM
        algorithm (more PCs = longer runtime).

    n_genes
        Determines the number of top SAM-weighted genes to use at each iteration
        of the SAM algorithm. If None, this number is chosen automatically
        based on the size of the dataset. This parameter primarily affects
        the runtime of the SAM algorithm (more genes = longer runtime). For
        extremely homogeneous datasets, decreasing `n_genes` may improve
        clustering resolution.

    inplace
        Set fields in `adata` if True. Otherwise, returns a copy.

    verbose
        If True, displays SAM log statements.

    Returns
    -------
    sam_obj if inplace is True or (sam_obj,AnnData) otherwise

    adata - AnnData
        `.var['weights']`
            SAM weights for each gene.
        `.var['spatial_dispersions']`
            Spatial dispersions for each gene (these are used to compute the
            SAM weights)
        `.uns['sam']`
            Dictionary of SAM-specific outputs, such as the parameters
            used for preprocessing ('preprocess_args') and running
            ('run_args') SAM.
        `.uns['neighbors']`
            A dictionary with key 'connectivities' containing the kNN adjacency
            matrix output by SAM. If built-in scanpy dimensionality reduction
            methods are to be used using the SAM-output AnnData, users
            should recompute the neighbors using `.obs['X_pca']` with
            `scanpy.pp.neighbors`.
        `.obsm['X_pca']`
            The principal components output by SAM.
        `.obsm['X_umap']`
            The UMAP projection output by SAM.
        `.layers['X_disp']`
            The expression matrix used for nearest-neighbor averaging.
        `.layers['X_knn_avg']`
            The nearest-neighbor-averaged expression data used for computing the
            spatial dispersions of genes.

    Example
    -------
    >>> import scanpy.external as sce
    >>> import scanpy as sc

    *** Running SAM ***

    Assuming we are given an AnnData object called `adata`, we can run the SAM
    algorithm as follows:

    >>> sam_obj = sce.tl.sam(adata,inplace=True)

    The input AnnData object should contain unstandardized, non-negative
    expression values. Preferably, the data should be log-normalized and no
    genes should be filtered out.

    Please see the documentation for a description of all available parameters.

    For more detailed tutorials, please visit the original Github repository:
    https://github.com/atarashansky/self-assembling-manifold/tree/master/tutorial

    *** Plotting ***

    To visualize the output, we can use:

    >>> sce.pl.sam(adata,projection='X_umap')

    `sce.pl.sam` accepts all keyword arguments used in the
    `matplotlib.pyplot.scatter` function.

    *** SAMGUI ***

    SAM comes with the SAMGUI module, a graphical-user interface written with
    `Plotly` and `ipythonwidgets` for interactively exploring and annotating
    the scRNAseq data and running SAM.

    Dependencies can be installed with Anaconda by following the instructions in
    the self-assembling-manifold Github README:
    https://github.com/atarashansky/self-assembling-manifold

    In a Jupyter notebook, execute the following to launch the interface:

    >>> from samalg.gui import SAMGUI
    >>> sam_gui = SAMGUI(sam_obj) # sam_obj is your SAM object
    >>> sam_gui.SamPlot

    This can also be enabled in Jupyer Lab by following the instructions in the
    self-assembling-manifold README.

    """

    try:
        from samalg import SAM
    except ImportError:
        raise ImportError(
            '\nplease install sam-algorithm: \n\n'
            '\tgit clone git://github.com/atarashansky/self-assembling-manifold.git\n'
            '\tcd self-assembling-manifold\n'
            '\tpip install .'
        )

    logg.info('Self-assembling manifold')

    s = SAM(counts=adata, inplace=inplace)

    logg.info('Running SAM')
    s.run(
        max_iter=max_iter,
        num_norm_avg=num_norm_avg,
        k=k,
        distance=distance,
        preprocessing=standardization,
        weight_PCs=weight_pcs,
        npcs=n_pcs,
        n_genes=n_genes,
        projection=projection,
        sparse_pca=sparse_pca,
        verbose=verbose,
    )

    s.adata.uns['sam'] = {}
    for attr in ['nnm', 'preprocess_args', 'run_args', 'ranked_genes']:
        s.adata.uns['sam'][attr] = s.adata.uns.pop(attr, None)

    return s if inplace else (s, s.adata)
