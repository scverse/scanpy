from .._settings import settings


def mnn_correct(*datas, var_index=None, var_subset=None, batch_key='batch', index_unique='-',
                batch_categories=None, k=20, sigma=1., cos_norm_in=True, cos_norm_out=True,
                svd_dim=None, var_adj=True, compute_angle=False, mnn_order=None, svd_mode='rsvd',
                do_concatenate=True, save_raw=False, n_jobs=None, **kwargs):
    """Correct batch effects by matching mutual nearest neighbors [Haghverdi18]_ [Kang18]_.

    This uses the implementation of `mnnpy
    <https://github.com/chriscainx/mnnpy>`__ [Kang18]_.

    Depending on `do_concatenate`, returns matrices or `AnnData` objects in the
    original order containing corrected expression values or a concatenated
    matrix or AnnData object.

    Be reminded that it is not advised to use the corrected data matrices for
    differential expression testing.

    More information and bug reports `here <https://github.com/chriscainx/mnnpy>`__.

    Parameters
    ----------
    datas : `numpy.ndarray` or :class:`~anndata.AnnData`
        Expression matrices or AnnData objects. Matrices should be shaped like
        n_obs * n_vars (n_cell * n_gene) and have consistent number of
        columns. AnnData objects should have same number of variables.
    var_index : `list` or `None`, optional (default: None)
        The index (list of str) of vars (genes). Necessary when using only a
        subset of vars to perform MNN correction, and should be supplied with
        var_subset. When datas are AnnData objects, var_index is ignored.
    var_subset : `list` or `None`, optional (default: None)
        The subset of vars (list of str) to be used when performing MNN
        correction. Typically, a list of highly variable genes (HVGs). When set
        to None, uses all vars.
    batch_key : `str`, optional (default: 'batch')
        The batch_key for AnnData.concatenate. Only valid when do_concatenate
        and supplying AnnData objects.
    index_unique : `str`, optional (default: '-')
        The index_unique for AnnData.concatenate. Only valid when do_concatenate
        and supplying AnnData objects.
    batch_categories : `list` or `None`, optional (default: None)
        The batch_categories for AnnData.concatenate. Only valid when
        do_concatenate and supplying AnnData objects.
    k : `int`, optional (default: 20)
        Number of mutual nearest neighbors.
    sigma : `float`, optional (default: 1)
        The bandwidth of the Gaussian smoothing kernel used to compute the
        correction vectors. Default is 1.
    cos_norm_in : `bool`, optional (default: True)
        Whether cosine normalization should be performed on the input data prior
        to calculating distances between cells.
    cos_norm_out : `bool`, optional (default: True)
        Whether cosine normalization should be performed prior to computing corrected expression values.
    svd_dim : `int` or `None`, optional (default: None)
        The number of dimensions to use for summarizing biological substructure
        within each batch. If None, biological components will not be removed
        from the correction vectors.
    var_adj : `bool`, optional (default: True)
        Whether to adjust variance of the correction vectors. Note this step
        takes most computing time.
    compute_angle : `bool`, optional (default: False)
        Whether to compute the angle between each cell’s correction vector and
        the biological subspace of the reference batch.
    mnn_order : `list` or `None`, optional (default: None)
        The order in which batches are to be corrected. When set to None, datas
        are corrected sequentially.
    svd_mode : `str`, optional (default: 'rsvd')
        One of 'svd', 'rsvd', and 'irlb'. 'svd' computes SVD using a
        non-randomized SVD-via-ID algorithm, while 'rsvd' uses a randomized
        version. 'irlb' performes truncated SVD by implicitly restarted Lanczos
        bidiagonalization (forked from https://github.com/airysen/irlbpy).
    do_concatenate : `bool`, optional (default: True)
        Whether to concatenate the corrected matrices or AnnData objects. Default is True.
    save_raw : `bool`, optional (default: False)
        Whether to save the original expression data in the .raw attribute of AnnData objects.
    n_jobs : `int` or `None`, optional (default: None)
        The number of jobs. When set to None, automatically uses the number of cores.
    kwargs : `dict` or `None`, optional (default: None)
        optional keyword arguments for irlb.

    Returns
    -------
    **datas** : :class:`~numpy.ndarray` or :class:`~anndata.AnnData`
        Corrected matrix/matrices or AnnData object/objects, depending on the
        input type and `do_concatenate`.
    **mnn_list** : ``List[pandas.DataFrame]``
        A list containing MNN pairing information as DataFrames in each iteration step.
    **angle_list** : ``List[Tuple[Optional[float], int]]`` or ``None``
        A list containing angles of each batch.
    """
    try:
        from mnnpy import mnn_correct as mnn_cor
        n_jobs = settings.n_jobs if n_jobs is None else n_jobs
        datas, mnn_list, angle_list = mnn_cor(
            *datas, var_index=var_index, var_subset=var_subset, batch_key=batch_key, index_unique=index_unique,
            batch_categories=batch_categories, k=k, sigma=sigma, cos_norm_in=cos_norm_in, cos_norm_out=cos_norm_out,
            svd_dim=svd_dim, var_adj=var_adj, compute_angle=compute_angle, mnn_order=mnn_order, svd_mode=svd_mode,
            do_concatenate=do_concatenate, save_raw=save_raw, n_jobs=n_jobs, **kwargs)
        return datas, mnn_list, angle_list
    except ImportError:
        raise ImportError(
            'Please install the package mnnpy '
            '(https://github.com/chriscainx/mnnpy). ')
