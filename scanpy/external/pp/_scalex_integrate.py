from typing import Optional

def scalex_integrate(
    adata, 
    batch_categories=None,
    batch_name='batch',
    outdir='scalex_output/',
    **kwargs
    ):
    """\
    Use scalex to integrate [Xiong22]_ to online integrate different experiments.

    SCALEX [Xiong22]_ is an online integration algorithm for integrating single-cell data from multiple experiments 
    and projecting new incoming datasets onto the existing cell-embedding space. This function uses the SCALEX from 
    ``scalex`` python package. The input to SCALEX can be filename or AnnData object, or list, or mixed of them. 
    SCALEX will do the preprocess if given the raw data, otherwise should set the ``processsed=True``
    
    Parameters
    ----------
    adata
        A path list of AnnData matrices to concatenate with. Each matrix is referred to as a 'batch'.
    batch_categories
        Categories for the batch annotation. By default, use increasing numbers.
    batch_name
        Use this annotation in obs as batches for training model. Default: 'batch'.
    outdir
        output AnnData will be saved as ``adata.h5ad`` in outdir if outdir is not None, default is None


    
    Returns
    -------
    The output folder contains:
    adata.h5ad
        The AnnData matrice after batch effects removal. The low-dimensional representation of the data is stored at adata.obsm['latent'].
    checkpoint
        model.pt contains the variables of the model and config.pt contains the parameters of the model.
    log.txt
        Records raw data information, filter conditions, model parameters etc.
    umap.pdf 
        UMAP plot for visualization.


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

    >>> 'X_scalex_umap' as well as 'X_umap' will in adata.obsm

    if adata has been preprocessed, please set the ``processed = True``
    >>> sce.pp.scalex_integrate(adata, processed=True)

    if use other meta information as ``batch``, e.g. sample_id
    >>> sce.pp.scalex_integrate(adata, batch_name = 'sample_id') 
    """
    try:
        import scalex
    except ImportError:
        raise ImportError("\nplease install scalex: \n\n\tpip install scalex")

    scalex.SCALEX(adata, batch_categories, batch_name=batch_name, outdir=outdir, **kwargs)

