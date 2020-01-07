"""\
Run the Self-Assembling Manifold algorithm
"""
from typing import Optional, Union

from anndata import AnnData

from ... import logging as logg


def sam(
    adata: AnnData,
    max_iter: int = 10,
    num_norm_avg: int = 50,
    k: int = 20,
    distance: str = 'correlation',
    standardization: Optional[str] = 'Normalizer',
    weight_pcs: bool = True,
    npcs: Optional[int]=None,
    n_genes: Optional[int]=None,
    projection: Optional[str]='umap',
    inplace: bool = True,
    verbose: bool = True

) -> Optional[AnnData]:
    """Self-Assembling Manifolds single-cell RNA sequencing analysis tool.

    SAM iteratively rescales the input gene expression matrix to emphasize
    genes that are spatially variable along the intrinsic manifold of the data.
    It outputs the gene weights, nearest neighbor matrix, and a 2D projection.

    The AnnData input should contain unstandardized, non-negative values.
    Preferably, the data should be log-normalized and no genes should be filtered out.


    Parameters
    ----------

    k - int, optional, default 20
        The number of nearest neighbors to identify for each cell.

    distance : string, optional, default 'correlation'
        The distance metric to use when identifying nearest neighbors.
        Can be any of the distance metrics supported by sklearn's 'pdist'.

    max_iter - int, optional, default 10
        The maximum number of iterations SAM will run.

    projection - str, optional, default 'umap'
        If 'tsne', generates a t-SNE embedding. If 'umap', generates a UMAP
        embedding. Otherwise, no embedding will be generated.

    standardization - str, optional, default 'Normalizer'
        If 'Normalizer', use sklearn.preprocessing.Normalizer, which
        normalizes expression data prior to PCA such that each cell has
        unit L2 norm. If 'StandardScaler', use
        sklearn.preprocessing.StandardScaler, which normalizes expression
        data prior to PCA such that each gene has zero mean and unit
        variance. Otherwise, do not normalize the expression data. We
        recommend using 'StandardScaler' for large datasets with many
        expected cell types and 'Normalizer' otherwise.

    num_norm_avg - int, optional, default 50
        The top 'num_norm_avg' dispersions are averaged to determine the
        normalization factor when calculating the weights. This prevents
        genes with large spatial dispersions from skewing the distribution
        of weights.

    weight_pcs - bool, optional, default True
        If True, scale the principal components by their eigenvalues. In
        datasets with many expected cell types, setting this to False might
        improve the resolution as these cell types might be encoded by low-
        variance principal components.

    npcs - int, optional, default None,
        Determines the number of top principal components selected at each
        iteration of the SAM algorithm. If None, this number is chosen
        automatically based on the size of the dataset. If weight_pcs is
        set to True, this parameter primarily affects the runtime of the SAM
        algorithm (more PCs = longer runtime).

    n_genes - int, optional, default None:
        Determines the number of top SAM-weighted genes to use at each iteration
        of the SAM algorithm. If None, this number is chosen automatically
        based on the size of the dataset. This parameter primarily affects
        the runtime of the SAM algorithm (more genes = longer runtime).
    
    inplace - bool, optional, default True:
        Set fields in `adata` if True. Otherwise, returns a copy.
    
    verbose - bool, optional, default True:
        If True, displays SAM log statements.

    Returns
    -------
    sam - SAM
        The SAM object

    adata - AnnData
        `.var['weights']`
            SAM weights for each gene.
        `.var['spatial_dispersions']`
            Spatial dispersions for each gene (these are used to compute the
            SAM weights)
        `.var['mask_genes']`
            If preprocessed with SAM, this boolean vector indicates which genes
            were filtered out (=False).
        `.uns['preprocess_args']`
            Dictionary of parameters used for preprocessing.
        `.uns['run_args']`
            Dictionary of parameters used for running SAM.
        `.uns['pca_obj']`
            The sklearn.decomposition.PCA object.
        `.uns['X_processed']`
            The standardized and SAM-weighted data fed into PCA.
        `.uns['neighbors']`
            A dictionary with key 'connectivities' containing the kNN adjacency
            matrix output by SAM. If built-in scanpy dimensionality reduction
            methods are to be used using the SAM-output AnnData, users
            should recompute the neighbors using `.obs['X_pca']` with
            `scanpy.pp.neighbors`.
        `.uns['ranked_genes']`
            Gene IDs ranked in descending order by their SAM weights.
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

    >>> sam,adata = sce.tl.SAM(adata,inplace=True)

    The input AnnData object should contain unstandardized, non-negative
    expression values. Preferably, the data should be log-normalized and no
    genes should be filtered out.

    Please see the documentation for a description of all available parameters.

    For more detailed tutorials, please visit the original Github repository:
    https://github.com/atarashansky/self-assembling-manifold/tree/master/tutorial

    *** Plotting ***

    To visualize the output, we can use the built-in `scatter` function (this
    assumes that `matplotlib` is installed.)

    >>> sam.scatter(projection = 'X_umap')

    `scatter` accepts all keyword arguments used in the
    `matplotlib.pyplot.scatter` function. Please visit the plotting tutorials
    for more information:

    https://github.com/atarashansky/self-assembling-manifold/tree/master/tutorial/SAM_Plotting

    *** SAMGUI ***

    SAM comes with the SAMGUI module, a graphical-user interface written with
    `Plotly` and `ipythonwidgets` for interactively exploring and annotating
    the scRNAseq data and running SAM.

    Dependencies can be installed with Anaconda by following the instructions in
    the self-assembling-manifold Github README:
    https://github.com/atarashansky/self-assembling-manifold

    In a Jupyter notebook, execute the following to launch the interface:

    >>> from SAMGUI import SAMGUI
    >>> sam_gui = SAMGUI(sam) # sam is your SAM object
    >>> sam_gui.SamPlot

    This can also be enabled in Jupyer Lab by following the instructions in the
    self-assembling-manifold README.

    """

    logg.info('Self-assembling manifold')

    try:
        from SAM import SAM
    except ImportError:
        raise ImportError(
            '\nplease install sam-algorithm: \n\n'
            '\tgit clone git://github.com/atarashansky/self-assembling-manifold.git\n'
            '\tcd self-assembling-manifold\n'
            '\tpip install .'
        )

    s = SAM(counts = adata, inplace = inplace)

    logg.info('Running SAM')
    s.run(max_iter=max_iter,num_norm_avg=num_norm_avg,k=k,
            distance=distance,preprocessing=standardization,
            weight_PCs=weight_pcs,npcs=npcs,n_genes=n_genes,
            projection = projection,verbose=verbose)

    return (s,adata) if inplace else (s,s.adata)

