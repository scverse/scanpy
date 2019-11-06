"""\
Run the Self-Assembling Manifold algorithm
"""
from typing import Optional, Union

from anndata import AnnData

from ... import logging as logg


def SAM(
    adata: AnnData,
    pp_with_sam: bool = True,
    pp_transform: Optional[str] = 'log',
    pp_min_expression: float = 1.0,
    pp_sum_norm: Union[str,float,None] = None,
    pp_thresh: float = 0.01,
    run_max_iter: int = 10,
    run_num_norm_avg: int = 50,
    run_k: int = 20,
    run_distance: str = 'correlation',
    run_standardization: Optional[str] = 'Normalizer',
    run_weight_PCs: bool = True,
    run_npcs: Optional[int]=None,
    run_n_genes: Optional[int]=None,
    run_projection: Optional[str]='umap'

) -> Optional[AnnData]:
    """Self-Assembling Manifolds single-cell RNA sequencing analysis tool.

    SAM iteratively rescales the input gene expression matrix to emphasize
    genes that are spatially variable along the intrinsic manifold of the data.
    It outputs the gene weights, nearest neighbor matrix, and a 2D projection.

    The AnnData input should contain unstandardized, non-negative values.

    Parameters with the 'pp_' prefix are data preprocessing parameters.
    Parameters with the 'run_' prefix are the SAM algorithm parameters.

    Parameters
    ----------

    pp_with_sam : bool, optional, default True
        If True, preprocess the data using built-in SAM functions. A copy
        of the input AnnData will be made in this case. If False, do not
        preprocess the data. SAM will be run in-place on the input AnnData.

    pp_sum_norm : str or float, optional, default None
        If a float, the total number of transcripts in each cell will be
        normalized to this value prior to normalization and filtering.
        Otherwise, nothing happens. If 'cell_median', each cell is
        normalized to have the median total read count per cell. If
        'gene_median', each gene is normalized to have the median total
        read count per gene.

    pp_transform : str, optional, default 'log'
        If 'log', log-transforms the expression data. If 'ftt', applies the
        Freeman-Tukey variance-stabilization transformation. If None, the data
        is not transformed.

    pp_min_expression : float, optional, default 1
        The threshold above which a gene is considered
        expressed. Gene expression values less than 'min_expression' are
        set to zero.

    pp_thresh : float, optional, default 0.2
        Keep genes expressed in greater than 'thresh'*100 % of cells and
        less than (1-'thresh')*100 % of cells, where a gene is considered
        expressed if its expression value exceeds 'min_expression'.

    run_k - int, optional, default 20
        The number of nearest neighbors to identify for each cell.

    run_distance : string, optional, default 'correlation'
        The distance metric to use when identifying nearest neighbors.
        Can be any of the distance metrics supported by sklearn's 'pdist'.

    run_max_iter - int, optional, default 10
        The maximum number of iterations SAM will run.

    run_projection - str, optional, default 'umap'
        If 'tsne', generates a t-SNE embedding. If 'umap', generates a UMAP
        embedding. Otherwise, no embedding will be generated.

    run_standardization - str, optional, default 'Normalizer'
        If 'Normalizer', use sklearn.preprocessing.Normalizer, which
        normalizes expression data prior to PCA such that each cell has
        unit L2 norm. If 'StandardScaler', use
        sklearn.preprocessing.StandardScaler, which normalizes expression
        data prior to PCA such that each gene has zero mean and unit
        variance. Otherwise, do not normalize the expression data. We
        recommend using 'StandardScaler' for large datasets with many
        expected cell types and 'Normalizer' otherwise.

    run_num_norm_avg - int, optional, default 50
        The top 'num_norm_avg' dispersions are averaged to determine the
        normalization factor when calculating the weights. This prevents
        genes with large spatial dispersions from skewing the distribution
        of weights.

    run_weight_PCs - bool, optional, default True
        If True, scale the principal components by their eigenvalues. In
        datasets with many expected cell types, setting this to False might
        improve the resolution as these cell types might be encoded by low-
        variance principal components.

    run_npcs - int, optional, default None,
        Determines the number of top principal components selected at each
        iteration of the SAM algorithm. If None, this number is chosen
        automatically based on the size of the dataset. If run_weight_PCs is
        set to True, this parameter primarily affects the runtime of the SAM
        algorithm (more PCs = longer runtime).

    run_n_genes - int, optional, default None:
        Determines the number of top SAM-weighted genes to use at each iteration
        of the SAM algorithm. If None, this number is chosen automatically
        based on the size of the dataset. This parameter primarily affects
        the runtime of the SAM algorithm (more genes = longer runtime).

    Returns
    -------
    sam - SAM
        The SAM object

    adata - AnnData
        The AnnData object. If pp_with_sam = True, this is simply a reference
        to the input AnnData object.

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

    >>> sam,adata = sce.tl.SAM(adata,pp_with_sam=True)

    The input AnnData object should contain unstandardized, non-negative
    expression values.

    Setting `pp_with_sam=True` preprocesses the data using a built-in function
    in SAM.

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

    if pp_with_sam:
        adata = adata.copy()
        sam = SAM(counts = adata)
        logg.info('Preprocessing data')
        sam.preprocess_data(sum_norm=pp_sum_norm, norm = pp_transform,
                            thresh = pp_thresh, min_expression=pp_min_expression)
    else:
        sam = SAM(counts = adata)

    logg.info('Running SAM')
    sam.run(max_iter=run_max_iter,num_norm_avg=run_num_norm_avg,k=run_k,
            distance=run_distance,preprocessing=run_standardization,
            weight_PCs=run_weight_PCs,npcs=run_npcs,n_genes=run_n_genes,
            projection = run_projection,verbose=False)

    return sam, sam.adata

