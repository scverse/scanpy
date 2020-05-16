import numpy as np
import pandas as pd
import scipy as sp

from typing import Optional, Sequence, Union
from anndata import AnnData

MIN_VERSION = "0.6.5"


def _prepare_dataset(
    adata: AnnData,
    layer: Optional[str] = None,
    subset_genes: Optional[Sequence[Union[int, str]]] = None,
    use_highly_variable_genes: bool = True,
    batch_key: Optional[str] = None,
    prepare_protein: bool = False,
):

    from scvi.dataset import AnnDatasetFromAnnData

    # check if observations are unnormalized using first 10
    # code from: https://github.com/theislab/dca/blob/89eee4ed01dd969b3d46e0c815382806fbfc2526/dca/io.py#L63-L69
    def _check_count_matrix(X):
        if len(adata) > 10:
            X_subset = X[:10]
        else:
            X_subset = X
        norm_error = 'Make sure that the dataset contains unnormalized count data.'
        if sp.sparse.issparse(X_subset):
            assert (X_subset.astype(int) != X_subset).nnz == 0, norm_error
        else:
            assert np.all(X_subset.astype(int) == X_subset), norm_error

    X = adata.X if layer is None else adata.layers[layer]
    _check_count_matrix(X)

    if subset_genes is not None:
        adata_subset = adata[:, subset_genes]
    elif use_highly_variable_genes and "highly_variable" in adata.var:
        adata_subset = adata[:, adata.var["highly_variable"]]
    else:
        adata_subset = adata

    adata_subset = adata_subset.copy()
    if layer is not None:
        adata_subset.X = adata_subset.layers[layer]
    if batch_key is not None:
        codes, uniques = pd.factorize(adata_subset.obs[batch_key])
        adata_subset.obs['_tmp_scvi_batch'] = codes

    if prepare_protein:
        _check_count_matrix(adata_subset.obsm["protein_expression"])
        mapping = {"protein_expression": "protein_names"}
    else:
        mapping = None
    dataset = AnnDatasetFromAnnData(
        adata_subset,
        batch_label='_tmp_scvi_batch',
        cell_measurements_col_mappings=mapping,
    )

    return dataset


def scvi(
    adata: AnnData,
    layer: Optional[str] = None,
    n_hidden: int = 128,
    n_latent: int = 10,
    n_layers: int = 1,
    dispersion: str = "gene",
    n_epochs: int = 400,
    lr: int = 1e-3,
    train_size: int = 1.0,
    batch_key: Optional[str] = None,
    use_highly_variable_genes: bool = True,
    subset_genes: Optional[Sequence[Union[int, str]]] = None,
    linear_decoder: bool = False,
    copy: bool = False,
    use_cuda: bool = True,
    return_posterior: bool = True,
    trainer_kwargs: dict = {},
    model_kwargs: dict = {},
) -> Optional[AnnData]:
    """\
    SCVI [Lopez18]_.

    Fits scVI model onto raw count data given an anndata object

    scVI uses stochastic optimization and deep neural networks to aggregate information
    across similar cells and genes and to approximate the distributions that underlie
    observed expression values, while accounting for batch effects and limited sensitivity.

    To use a linear-decoded Variational AutoEncoder model (implementation of [Svensson20]_.),
    set linear_decoded = True. Compared to standard VAE, this model is less powerful, but can
    be used to inspect which genes contribute to variation in the dataset. It may also be used
    for all scVI tasks, like differential expression, batch correction, imputation, etc.
    However, batch correction may be less powerful as it assumes a linear model.

    .. note::
        This is a limited application of totalVI.
        Full API and more information `here <https://scvi.readthedocs.io/en/stable/>`__.

    Parameters
    ----------
    adata
        An anndata file with `X` attribute of unnormalized count data
    layer
        Layer of `adata` to use as expression values.
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dispersion
        One of the following

        * `'gene'` - dispersion parameter of NB is constant per gene across cells
        * `'gene-batch'` - dispersion can differ between different batches
        * `'gene-label'` - dispersion can differ between different labels
        * `'gene-cell'` - dispersion can differ for every gene in every cell
    n_epochs
        Number of epochs to train
    lr
        Learning rate
    train_size
        The train size, either a float between 0 and 1 or an integer for the number of training samples to use
    batch_key
        Column name in anndata.obs for batches.
        If None, no batch correction is performed
        If not None, batch correction is performed per batch category
    use_highly_variable_genes
        If true, uses only the genes in `anndata.var["highly_variable"]`
    subset_genes
        Optional list of indices or gene names to subset anndata.
        If not None, use_highly_variable_genes is ignored
    linear_decoder
        If true, uses LDVAE model, which is an implementation of [Svensson20]_.
    copy
        If true, a copy of anndata is returned
    return_posterior
        If true, posterior object is returned
    use_cuda
        If true, uses cuda
    trainer_kwargs
        Extra arguments for UnsupervisedTrainer
    model_kwargs
        Extra arguments for VAE or LDVAE model

    Returns
    -------
    If `copy` is true, anndata is returned.
    If `return_posterior` is true, the posterior object is returned
    If both `copy` and `return_posterior` are true,
    a tuple of anndata and the posterior are returned in that order.

    By default the following items are stored in the input adata:

    `adata.obsm['X_scvi']` stores the latent representations
    `adata.obsm['X_scvi_denoised']` stores the normalized mean of the negative binomial
    `adata.obsm['X_scvi_sample_rate']` stores the mean of the negative binomial

    If linear_decoder is true:
    `adata.uns['ldvae_loadings']` stores the per-gene weights in the linear decoder as a
    genes by n_latent matrix.

    Examples
    --------
    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> adata = sc.datasets.pbmc3k()
    >>> adata.var_names_make_unique()
    >>> adata.layers['counts'] = adata.X
    >>> sc.pp.normalize_per_cell(adata)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata)
    >>> posterior = sce.pp.scvi(adata, layer='counts', return_posterior=True)
    >>> sc.pp.neighbors(adata, use_rep="X_scvi")
    >>> sc.tl.umap(adata)
    """

    try:
        import scvi as sc
        from scvi.models import VAE, LDVAE
        from scvi.inference import UnsupervisedTrainer
    except ImportError:
        raise ImportError(
            "Please install scvi package from https://github.com/YosefLab/scVI"
        )

    if sc.__version__ < MIN_VERSION:
        raise ValueError("Please upgrade scvi via `pip install --upgrade scvi`")

    dataset = _prepare_dataset(
        adata,
        layer=layer,
        subset_genes=subset_genes,
        use_highly_variable_genes=use_highly_variable_genes,
        batch_key=batch_key,
    )

    if linear_decoder:
        vae = LDVAE(
            n_input=dataset.nb_genes,
            n_batch=dataset.n_batches,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers,
            dispersion=dispersion,
            **model_kwargs,
        )

    else:
        vae = VAE(
            dataset.nb_genes,
            n_batch=dataset.n_batches,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dispersion=dispersion,
            **model_kwargs,
        )

    trainer = UnsupervisedTrainer(
        model=vae,
        gene_dataset=dataset,
        use_cuda=use_cuda,
        train_size=train_size,
        **trainer_kwargs,
    )

    trainer.train(n_epochs=n_epochs, lr=lr)

    full = trainer.create_posterior()
    latent, _, _ = full.get_latent()

    if copy:
        adata = adata.copy()

    adata.obsm['X_scvi'] = latent
    adata.obsm['X_scvi_denoised'] = pd.DataFrame(
        full.get_sample_scale(), columns=dataset.gene_names, index=adata.obs_names,
    )
    adata.obsm['X_scvi_sample_rate'] = pd.DataFrame(
        full.imputation(), columns=dataset.gene_names, index=adata.obs_names,
    )

    if linear_decoder:
        loadings = vae.get_loadings()
        df = pd.DataFrame(loadings, index=dataset.gene_names)
        adata.uns['ldvae_loadings'] = df

    if copy and return_posterior:
        return adata, full
    elif copy:
        return adata
    elif return_posterior:
        return full


def totalvi(
    adata: AnnData,
    layer: Optional[str] = None,
    n_latent: int = 20,
    n_epochs: int = 500,
    n_hidden: int = 256,
    lr: int = 4e-3,
    train_size: int = 0.9,
    early_stopping: bool = True,
    batch_key: Optional[str] = None,
    use_highly_variable_genes: bool = True,
    subset_genes: Optional[Sequence[Union[int, str]]] = None,
    copy: bool = False,
    use_cuda: bool = True,
    return_posterior: bool = True,
    trainer_kwargs: dict = {},
    model_kwargs: dict = {},
) -> Optional[AnnData]:
    """\
    TOTALVI [Gayoso20]_.

    Fits totalVI model onto raw CITE-seq count data given an anndata object

    totalVI, like scVI, is based on a hierarchical Bayesian model and uses the framework of
    variational autoencoders to model CITE-seq data (protein + RNA). Special attention is given
    to the protein data, which is biased by background due to ambient antibodies/non-specific binding.

    .. note::
        This is a limited application of totalVI.
        Full API and more information `here <https://scvi.readthedocs.io/en/stable/>`__.

    Parameters
    ----------
    adata
        An anndata file with `X` attribute of unnormalized count data and
        `protein_expression` key in `.obsm` of cells by proteins count matrix and
        `protein_names` key in `.uns` that is 1d array of str with protein names
        Both `protein_expression` and `protein_names` are expected to be (:class:`~numpy.ndarray`)
    layer
        Layer of `adata` to use as expression values.
    n_latent
        Dimensionality of the latent space
    n_epochs
        Number of epochs to train
    n_hidden
        Number of nodes per hidden layer
    lr
        Learning rate
    train_size
        The train size, either a float between 0 and 1 or an integer for the number of training samples to use
    early_stopping
        Whether to use early stopping during training. If `False`, early stopping can still be used by
        adding custom `early_stopping_kwargs` to `trainer_kwargs`
    batch_key
        Column name in anndata.obs for batches.
        If None, no batch correction is performed
        If not None, batch correction is performed per batch category
    use_highly_variable_genes
        If true, uses only the genes in `anndata.var["highly_variable"]`
    subset_genes
        Optional list of indices or gene names to subset anndata.
        If not None, use_highly_variable_genes is ignored
    copy
        If true, a copy of anndata is returned
    return_posterior
        If true, posterior object is returned
    use_cuda
        If true, uses cuda
    trainer_kwargs
        Extra arguments for TotalTrainer
    model_kwargs
        Extra arguments for TOTALVI model

    Returns
    -------
    If `copy` is true, anndata is returned.
    If `return_posterior` is true, the posterior object is returned
    If both `copy` and `return_posterior` are true,
    a tuple of anndata and the posterior are returned in that order.

    `adata.obsm['X_totalvi']` stores the latent representations
    `adata.obsm['X_totalvi_denoised_genes']` stores the normalized mean of the negative binomial
    `adata.obsm['X_totalvi_denoised_proteins']` stores the background corrected protein values
    `adata.obsm['X_totalvi_protein_prob_foreground']` stores the probability a protein measurement is
    from the "foreground" (i.e., not due to background antibodies)

    """

    try:
        import scvi as sc
        from scvi.models import TOTALVI
        from scvi.inference import TotalPosterior, TotalTrainer
    except ImportError:
        raise ImportError(
            "Please install scvi package from https://github.com/YosefLab/scVI"
        )
    if sc.__version__ < MIN_VERSION:
        raise ValueError("Please upgrade scvi via `pip install --upgrade scvi`")

    dataset = _prepare_dataset(
        adata,
        layer=layer,
        subset_genes=subset_genes,
        use_highly_variable_genes=use_highly_variable_genes,
        batch_key=batch_key,
        prepare_protein=True,
    )

    vae = TOTALVI(
        dataset.nb_genes,
        dataset.protein_expression.shape[1],
        n_batch=dataset.n_batches,
        n_hidden=n_hidden,
        n_latent=n_latent,
        **model_kwargs,
    )

    if early_stopping is True:
        trainer_kwargs.update({'early_stopping_kwargs': 'auto'})
        frequency = 1
    else:
        trainer_kwargs.update({'early_stopping_kwargs': None})
        frequency = None

    trainer = TotalTrainer(
        model=vae,
        dataset=dataset,
        use_cuda=use_cuda,
        train_size=train_size,
        test_size=1.0 - train_size,
        frequency=frequency,
        **trainer_kwargs,
    )

    trainer.train(n_epochs=n_epochs, lr=lr)

    full = trainer.create_posterior(type_class=TotalPosterior)
    latent = full.get_latent()[0]

    if copy:
        adata = adata.copy()

    adata.obsm['X_totalvi'] = latent
    denoised_genes, denoised_proteins = full.get_normalized_denoised_expression(
        n_samples=10, give_mean=True
    )
    adata.obsm['X_totalvi_denoised_genes'] = pd.DataFrame(
        denoised_genes, index=adata.obs_names, columns=dataset.gene_names
    )
    adata.obsm['X_totalvi_denoised_proteins'] = pd.DataFrame(
        denoised_proteins, index=adata.obs_names, columns=dataset.protein_names
    )
    adata.obsm['X_totalvi_protein_prob_foreground'] = pd.DataFrame(
        1.0 - full.get_sample_mixing(n_samples=10, give_mean=True),
        index=adata.obs_names,
        columns=dataset.protein_names,
    )

    if copy and return_posterior:
        return adata, full
    elif copy:
        return adata
    elif return_posterior:
        return full
