import numpy as np
import pandas as pd
import scipy as sp

from typing import Optional, Sequence, Union
from anndata import AnnData

MIN_VERSION = "0.6.5"


def scvi(
    adata: AnnData,
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
        More information and bug reports `here <https://github.com/YosefLab/scVI>`__.

    Parameters
    ----------
    adata
        An anndata file with `X` attribute of unnormalized count data
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
        If true, uses only the genes in anndata.var["highly_variable"]
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

    `adata.obsm['X_scvi']` stores the latent representations
    `adata.obsm['X_scvi_denoised']` stores the normalized mean of the negative binomial
    `adata.obsm['X_scvi_sample_rate']` stores the mean of the negative binomial
    
    If linear_decoder is true:
    `adata.uns['ldvae_loadings']` stores the per-gene weights in the linear decoder as a
    genes by n_latent matrix.

    """

    try:
        from scvi.models import VAE, LDVAE
        from scvi.inference import UnsupervisedTrainer
        from scvi.dataset import AnnDatasetFromAnnData
    except ImportError:
        raise ImportError(
            "Please install scvi package from https://github.com/YosefLab/scVI"
        )

    # check if observations are unnormalized using first 10
    # code from: https://github.com/theislab/dca/blob/89eee4ed01dd969b3d46e0c815382806fbfc2526/dca/io.py#L63-L69
    if len(adata) > 10:
        X_subset = adata.X[:10]
    else:
        X_subset = adata.X
    norm_error = (
        'Make sure that the dataset (adata.X) contains unnormalized count data.'
    )
    if sp.sparse.issparse(X_subset):
        assert (X_subset.astype(int) != X_subset).nnz == 0, norm_error
    else:
        assert np.all(X_subset.astype(int) == X_subset), norm_error

    if subset_genes is not None:
        adata_subset = adata[:, subset_genes]
    elif use_highly_variable_genes and "highly_variable" in adata.var:
        adata_subset = adata[:, adata.var["highly_variable"]]
    else:
        adata_subset = adata

    if batch_key is not None:
        codes, uniques = pd.factorize(adata_subset.obs[batch_key])
        adata_subset.obs['_tmp_scvi_batch'] = codes
        n_batches = len(uniques)
    else:
        n_batches = 0

    dataset = AnnDatasetFromAnnData(adata_subset.copy(), batch_label='_tmp_scvi_batch')

    if linear_decoder:
        vae = LDVAE(
            n_input=dataset.nb_genes,
            n_batch=n_batches,
            n_labels=dataset.n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers,
            dispersion=dispersion,
            **model_kwargs,
        )

    else:
        vae = VAE(
            dataset.nb_genes,
            n_batch=n_batches,
            n_labels=dataset.n_labels,
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

    full = trainer.create_posterior(
        trainer.model, dataset, indices=np.arange(len(dataset))
    )
    latent, batch_indices, labels = full.sequential().get_latent()

    if copy:
        adata = adata.copy()

    adata.obsm['X_scvi'] = latent
    adata.obsm['X_scvi_denoised'] = full.sequential().get_sample_scale()
    adata.obsm['X_scvi_sample_rate'] = full.sequential().imputation()

    if linear_decoder:
        loadings = vae.get_loadings()
        df = pd.DataFrame(loadings, index=adata_subset.var_names)
        adata.uns['ldvae_loadings'] = df

    if copy and return_posterior:
        return adata, full
    elif copy:
        return adata
    elif return_posterior:
        return full
