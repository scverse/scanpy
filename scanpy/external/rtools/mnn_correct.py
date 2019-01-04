# This is not actually used as we have a wrapper for a Python package in the meanwhile...
# It might just serve as a template for future integrations...

def mnn_concatenate(*adatas, geneset=None, k=20, sigma=1, n_jobs=None, **kwds):
    """Merge AnnData objects and correct batch effects using the MNN method.

    Batch effect correction by matching mutual nearest neighbors [Haghverdi18]_
    has been implemented as a function 'mnnCorrect' in the R package
    `scran <https://bioconductor.org/packages/release/bioc/html/scran.html>`__
    This function provides a wrapper to use the mnnCorrect function when
    concatenating Anndata objects by using the Python-R interface `rpy2
    <https://pypi.org/project/rpy2/>`__.

    Parameters
    ----------
    adatas : :class:`~anndata.AnnData`
        AnnData matrices to concatenate with. Each dataset should generally be
        log-transformed, e.g., log-counts. Datasets should have the same number
        of genes, or at lease have all the genes in geneset.
    geneset : `list`, optional (default: `None`)
        A list specifying the genes with which distances between cells are
        calculated in mnnCorrect, typically the highly variable genes.
        All genes are used if no geneset provided. See the `scran manual
        <https://bioconductor.org/packages/release/bioc/html/scran.html>`__ for
        details.
    k : `int`, ptional (default: 20)
        See the `scran manual <https://bioconductor.org/packages/release/bioc/html/scran.html>`__
        for details.
    sigma : `int`, ptional (default: 20)
        See the `scran manual <https://bioconductor.org/packages/release/bioc/html/scran.html>`__
        for details.
    n_jobs : `int` or `None` (default: `sc.settings.n_jobs`)
        Number of jobs.
    kwds :
        Keyword arguments passed to Anndata.concatenate

    Returns
    -------
    An :class:`~anndata.AnnData` object with MNN corrected data matrix X.

    Example
    -------
    >>> adata1
    AnnData object with n_obs × n_vars = 223 × 33694
        obs: 'n_genes', 'percent_mito', 'n_counts', 'Sample', 'Donor', 'Tissue'
        var: 'gene_ids', 'n_cells'
    >>> adata2
    AnnData object with n_obs × n_vars = 1457 × 33694
        obs: 'n_genes', 'percent_mito', 'n_counts', 'Sample', 'Donor', 'Tissue'
        var: 'gene_ids', 'n_cells'
    >>> adata3 = sc.pp.mnnconcatenate(adata2, adata1, geneset = hvgs)
    """
    from rpy2.robjects.packages import importr
    from rpy2.robjects import numpy2ri

    adata = AnnData.concatenate(*adatas, **kwds)
    if geneset is None:
        datamats = tuple([adata.X.T for adata in adatas])
    else:
        datamats = tuple([adata[:, geneset].X.T for adata in adatas])
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs
    numpy2ri.activate()
    rbase      = importr('base')
    rscran     = importr('scran')
    bpparam    = importr('BiocParallel').MulticoreParam(
        workers = n_jobs) if n_jobs > 1 else importr('BiocParallel').SerialParam()
    mnn_result = rscran.mnnCorrect(*datamats, k=k, sigma=sigma, BPPARAM = bpparam)
    corrected  = np.asarray(rbase.do_call(rbase.cbind, mnn_result[0])).T
    numpy2ri.deactivate()
    if geneset is None:
        adata = adata[:, geneset]
    adata.X = corrected
    return adata

