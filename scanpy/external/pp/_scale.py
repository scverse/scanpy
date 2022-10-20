from typing import Optional, Union, List
from anndata import AnnData


def scale(
    data_list: Union[str, List],
    n_centroids: int = 30,
    outdir: bool = None,
    verbose: bool = False,
    pretrain: str = None,
    lr: float = 0.0002,
    batch_size: int = 64,
    gpu: int = 0,
    seed: int = 18,
    encode_dim: List = [1024, 128],
    decode_dim: List = [],
    latent: int = 10,
    min_peaks: int = 100,
    min_cells: Union[float, int] = 3,
    n_feature: int = 50000,
    log_transform: bool = False,
    max_iter: int = 30000,
    weight_decay: float = 5e-4,
    impute: bool = False,
    binary: bool = False,
    embed: str = 'UMAP',
    reference: str = 'cell_type',
    cluster_method: str = 'leiden',
) -> AnnData:
    """
    Use SCALE [Xiong19]_ to cluster and impute on scATAC-seq

    SCALE [Xiong19]_ is a deeplearning tool combining GMM with VAE for single-cell ATAC-seq analysis (visualization, clustering, imputation, downstream analysis for celltype-specific TFs

    Parameters
    ----------
    data_list
        A path of input data, could be AnnData, or file in format of h5ad, txt, csv or dir containing mtx file
    n_centroids
        Initial n_centroids, default is 30
    outdir
        output of dir where results will be saved, if None, will return an AnnData object
    verbose
        Verbosity, default is False
    pretrain
        Use the pretrained model to project on new data or repeat results
    lr
        learning rate for training model
    batch_size
        batch_size of one iteration for training model
    gpu
        Use id of gpu device
    seed
        Random seed
    encode_dim
        encode architecture
    decode_dim
        decode architecture
    latent
        dimensions of latent
    min_peaks
        min peaks for filtering cells
    min_cells
        min cells for filtering peaks, will be ratio if small than 1, default is 3
    n_feature
        number of the most highly variable peaks will be used for training model
    log_transform
        log transform the data, recommended for scRNA-seq data
    max_iter
        Max iteration for training model
    weight_decay
        weight decay for training model
    impute
        output the imputed data in adata if true
    binary
        Change the imputed data to binary format
    embed
        Embedding method, UMAP or tSNE, default is UMAP
    reference
        reference annotations for evaluation, default is cell_type
    cluster_method
        cluster method, default is Leiden
    )->AnnData:


    Returns
    -------
    adata
        An AnnData object


    Example
    -------
    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> adata = sc.read('Your/scATAC/Data')
    >>> adata = sce.pp.scale(adata)

    if want imputed data

    >>> adata = sce.pp.scale(adata, impute=True)

    imputed data will be saved in adata.obsm['impute']
    binary version of imputed data will be saved in adata.obsm['binary']
    """

    try:
        import scale
    except ImportError:
        raise ImportError("\nplease install scale: \n\n\tpip install scale")

    return scale.SCALE_function(
        data_list,
        n_centroids=n_centroids,
        outdir=outdir,
        verbose=verbose,
        pretrain=pretrain,
        lr=lr,
        batch_size=batch_size,
        gpu=gpu,
        seed=seed,
        encode_dim=encode_dim,
        decode_dim=decode_dim,
        latent=latent,
        min_peaks=min_peaks,
        min_cells=min_cells,
        n_feature=n_feature,
        log_transform=log_transform,
        max_iter=max_iter,
        weight_decay=weight_decay,
        impute=impute,
        binary=binary,
        embed=embed,
        reference=reference,
        cluster_method=cluster_method,
    )
