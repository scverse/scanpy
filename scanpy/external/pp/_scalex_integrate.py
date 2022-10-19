from typing import Optional, Union, List
from anndata import AnnData


def scalex_integrate(
    data_list: Union[str, AnnData, List] = None,
    batch_categories: List = None,
    profile: str = 'RNA',
    batch_name: str = 'batch',
    min_features: int = 600,
    min_cells: int = 3,
    target_sum: int = None,
    n_top_features: int = None,
    join: str = 'inner',
    batch_key: str = 'batch',
    processed: bool = False,
    batch_size: int = 64,
    lr: float = 2e-4,
    max_iteration: int = 30000,
    seed: int = 124,
    gpu: int = 0,
    outdir: str = None,
    projection: str = None,
    repeat: bool = False,
    impute: str = None,
    chunk_size: int = 20000,
    ignore_umap: bool = False,
    verbose: bool = False,
    assess: bool = False,
    show: bool = True,
    eval: bool = False,
) -> AnnData:
    """
    Use scalex [Xiong22]_ to online integrate different experiments.

    SCALEX [Xiong22]_ is an algorithm for online integrating single-cell
    data from multiple experiments. This function uses the ``scalex.SCALEX``, to integrate single-cell data
    stored in an AnnData object, or list of file paths. SCALEX prefer to take the raw input data without any preprocessing.
    If the data has already been preprocessed, should set `processed=True`.
    Online single-cell data integration through projecting heterogeneous datasets into a common cell-embedding space

    Parameters
    ----------
    data_list
        A path list of AnnData matrices to concatenate with. Each matrix is referred to as a 'batch'.
    batch_categories
        Categories for the batch annotation. By default, use increasing numbers.
    profile
        Specify the single-cell profile, RNA or ATAC. Default: RNA.
    batch_name
        Use this annotation in obs as batches for training model. Default: 'batch'.
    min_features
        Filtered out cells that are detected in less than min_features. Default: 600.
    min_cells
        Filtered out genes that are detected in less than min_cells. Default: 3.
    n_top_features
        Number of highly-variable genes to keep. Default: 2000.
    join
        Use intersection ('inner') or union ('outer') of variables of different batches.
    batch_key
        Add the batch annotation to obs using this key. By default, batch_key='batch'.
    batch_size
        Number of samples per batch to load. Default: 64.
    lr
        Learning rate. Default: 2e-4.
    max_iteration
        Max iterations for training. Training one batch_size samples is one iteration. Default: 30000.
    seed
        Random seed for torch and numpy. Default: 124.
    gpu
        Index of GPU to use if GPU is available. Default: 0.
    outdir
        Output directory. Default: 'output/'.
    projection
        Use for new dataset projection. Input the folder containing the pre-trained model. If None, don't do projection. Default: None.
    repeat
        Use with projection. If False, concatenate the reference and projection datasets for downstream analysis. If True, only use projection datasets. Default: False.
    impute
        If True, calculate the imputed gene expression and store it at adata.layers['impute']. Default: False.
    chunk_size
        Number of samples from the same batch to transform. Default: 20000.
    ignore_umap
        If True, do not perform UMAP for visualization and leiden for clustering. Default: False.
    verbose
        Verbosity, True or False. Default: False.
    assess
        If True, calculate the entropy_batch_mixing score and silhouette score to evaluate integration results. Default: False.



    Returns
    -------
    adata
        An AnnData object


    Example
    -------
    >>> import scanpy as sc
    >>> import scapy.external as sce
    >>> adata = sc.datasets.pbmc3k()

    We now arbitrarily assign a batch metadata variable to each cell
    for the sake of example, but during real usage there would already
    be a column in ``adata.obs`` giving the experiment each cell came
    from.

    >>> adata.obs['batch] = 1350*['a'] + 1350*['b']
    >>> sce.pp.scalex_integrate(adata)

    'X_scalex_umap' as well as 'X_umap' will in adata.obsm

    SCALEX will do the preprocssing internally, like sc.pp.filter_cells, sc.pp.filter_genes, sc.pp.normalize_total, sc.pp.log1p, sc.pp.highly_variable_genes
    if adata has been preprocessed, please set the ``processed = True``

    >>> sce.pp.scalex_integrate(adata, processed=True)

    if use other meta information as ``batch``, e.g. sample_id

    >>> sce.pp.scalex_integrate(adata, batch_name = 'sample_id')
    """

    try:
        import scalex
    except ImportError:
        raise ImportError("\nplease install scalex: \n\n\tpip install scalex")

    return scalex.SCALEX(
        data_list=data_list,
        batch_categories=batch_categories,
        profile=profile,
        batch_name=batch_name,
        min_features=min_features,
        min_cells=min_cells,
        target_sum=target_sum,
        n_top_features=n_top_features,
        join=join,
        batch_key=batch_key,
        processed=processed,
        batch_size=batch_size,
        lr=lr,
        max_iteration=max_iteration,
        seed=seed,
        gpu=gpu,
        outdir=outdir,
        projection=projection,
        repeat=repeat,
        impute=impute,
        chunk_size=chunk_size,
        ignore_umap=ignore_umap,
        verbose=verbose,
        assess=assess,
        show=show,
        eval=eval,
    )
