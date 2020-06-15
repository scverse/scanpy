import collections
import numpy as np
import scipy.sparse as sp

from typing import Optional, Sequence
from anndata import AnnData

from ...preprocessing import filter_genes
from ..rtools import _check_rpy2, _check_anndata2ri, _is_installed, _py2r, _r2py, _set_seed


def sctransform(
    adata: AnnData,
    regress_out: Sequence = ('log_umi',),
    batch_key: Optional[str] = None,
    n_top_genes: int = 3000,
    min_cells: int = 5,
    store_residuals: bool = False,
    verbose: bool = True,
    inplace: bool = True,
    seed: int = 0,
) -> Optional[AnnData]:
    """\
    Normalization and variance stabilization of scRNA-seq data using regularized
    negative binomial regression [Hafemeister19]_.

    sctransform uses Pearson residuals from “regularized negative binomial regression,” to
    correct for the sequencing depth. After regressing out total number of UMIs (and other
    variables if given) it ranks the genes based on their residual variances and therefore
    also acts as a HVG selection method.

    This function replaces `sc.pp.normalize_total` and `sc.pp.highly_variable_genes` and requires
    raw counts in `adata.X`.

    .. note::
        More information and bug reports `here <https://github.com/ChristophH/sctransform>`__.

    Parameters
    ----------
    adata
        An anndata file with `X` attribute of unnormalized count data
    regress_out
        Variables to regress out. Default is logarithmized total UMIs which is implicitly
        calculated by sctransform. Other obs keys can also be used.
    batch_key
        If specified, HVGs are ranked after batch_key is regressed out. This avoids the
        selection of batch-specific genes and acts as a lightweight batch correction method.
        Note that corrected counts are not batch-corrected but only depth-normalized.
    n_top_genes
        Total number of highly variable genes selected.
    min_cells
        Only use genes that have been detected in at least this many cells; default is 5.
    store_residuals
        Store Pearson residuals in adata.layers['sct_residuals']. These values represent
        batch corrected and depth-normalized gene expression values. Due to potential
        high memory use for big matrices, they are not stored by default.
    verbose
        Show progress bar during normalization.
    inplace
        Save HVGs and corrected UMIs inplace. Default is True.
    seed
        Random seed for R RNG. Default is 0.

    Returns
    -------
    If `inplace` is False, anndata is returned.
    If `store_residuals` is True, residuals are stored in adata.layers['sct_residuals'].

    `adata.layers['sct_corrected']` stores normalized representation of gene expression.
    `adata.var['highly_variable']` stores highly variable genes.
    `adata.var['highly_variable_sct_residual_var']` stores the residual variances that
    are also used for ranking genes by variability.

    """

    _check_rpy2()
    _check_anndata2ri()

    import rpy2
    from rpy2.robjects import r
    from rpy2.robjects.packages import importr
    import anndata2ri

    _is_installed('sctransform')
    _set_seed(seed)

    # check if observations are unnormalized using first 10
    if len(adata) > 10:
        X_subset = adata.X[:10]
    else:
        X_subset = adata.X
    err = 'Make sure that adata.X contains unnormalized count data'
    if sp.issparse(X_subset):
        assert (X_subset.astype(int) != X_subset).nnz == 0, err
    else:
        assert np.all(X_subset.astype(int) == X_subset), err

    assert regress_out, 'regress_out cannot be emtpy'

    if not inplace:
        adata = adata.copy()

    filter_genes(adata, min_cells=min_cells)

    mat = adata.X.T
    if sp.issparse(mat):
        mat.sort_indices()
    mat = _py2r(mat)

    set_colnames = r('`colnames<-`')
    set_rownames = r('`rownames<-`')

    mat = set_colnames(mat, adata.obs_names.values.tolist())
    mat = set_rownames(mat, adata.var_names.values.tolist())

    assert isinstance(
        regress_out, collections.abc.Sequence
    ), 'regress_out is not a sequence'
    obs_keys = [x for x in regress_out if x != 'log_umi']
    regress_out = _py2r(np.array(regress_out))

    if batch_key is not None:
        obs_keys += [batch_key]
    else:
        batch_key = rpy2.rinterface.NULL

    if obs_keys:
        assert np.all(
            np.isin(obs_keys, adata.obs.columns)
        ), 'Some regress_out or batch_key values are not found in adata.obs'
        cell_attr = adata.obs[obs_keys]
        cell_attr = _py2r(cell_attr)
    else:
        cell_attr = rpy2.rinterface.NULL

    sct = importr('sctransform')
    residual_type = 'pearson' if store_residuals else 'none'

    vst_out = sct.vst(
        mat,
        cell_attr=cell_attr,
        batch_var=batch_key,
        latent_var=regress_out,
        residual_type=residual_type,
        return_cell_attr=True,
        min_cells=min_cells,
        show_progress=verbose,
    )

    res_var = _r2py(sct.get_residual_var(vst_out, mat))
    corrected = _r2py(sct.correct_counts(vst_out, mat)).T

    adata.layers['sct_corrected'] = corrected
    adata.var['highly_variable_sct_residual_var'] = res_var

    if store_residuals:
        adata.layers['sct_residuals'] = _r2py(vst_out.rx2('y')).T

    top_genes = (
        adata.var['highly_variable_sct_residual_var']
        .sort_values(ascending=False)[:n_top_genes]
        .index.tolist()
    )
    adata.var['highly_variable'] = adata.var_names.isin(top_genes)

    if not inplace:
        return adata
