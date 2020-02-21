import numpy as np

from typing import Optional, List
from anndata import AnnData


def scvi(
    adata: AnnData,
    n_hidden: int = 128,
    n_latent: int = 30,
    n_layers: int = 2,
    dispersion: str = "gene",
    n_epochs: int = 1,
    train_size: int = 1.0,
    use_cuda: bool = True,
    metrics_to_monitor: List = None,
    metrics_monitor_frequency: int = None,
    return_posterior: bool = False,
    return_info: bool = False,
    copy: bool = False,
    batch_label: str = 'batch_indices',
    ctype_label: str = 'cell_types',
    class_label: str = 'labels',
) -> Optional[AnnData]:
    """\
    SCVI [Lopez18]_.

    Fits SCVI model onto raw count data given in the anndata object

    .. note::
        More information and bug reports `here <https://github.com/YosefLab/scVI>`__.

    Parameters
    ----------
    adata
        An anndata file with `X` attribute representing raw counts.
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dispersion
        One of the following
        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    n_epochs: int = 1,
    train_size: int = 1.0,
    use_cuda: bool = True,
    metrics_to_monitor: List = None,
    metrics_monitor_frequency: int = None,
    return_posterior: bool = False,
    return_info: bool = False,
    copy: bool = False,
    batch_label: str = 'batch_indices',
    ctype_label: str = 'cell_types',
    class_label: str = 'labels',
    """

    try:
        from scvi.models import VAE
        from scvi.inference import UnsupervisedTrainer
        from scvi.dataset import AnnDatasetFromAnnData
    except ImportError:
        raise ImportError("Please install scvi package via")

    dataset = AnnDatasetFromAnnData(adata, use_raw=True)

    vae = VAE(
        dataset.nb_genes,
        n_batch=dataset.n_batches,
        n_labels=dataset.n_labels,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dispersion=dispersion,
    )

    trainer = UnsupervisedTrainer(
        model=vae,
        gene_dataset=dataset,
        use_cuda=use_cuda,
        metrics_to_monitor=metrics_to_monitor,
        frequency=metrics_monitor_frequency,
        train_size=train_size,
    )

    trainer.train(n_epochs=n_epochs)

    full = trainer.create_posterior(
        trainer.model, dataset, indices=np.arange(len(dataset))
    )
    latent, batch_indices, labels = full.sequential().get_latent()

    if copy:
        adata = adata.copy()

    adata.obsm['X_scvi'] = latent
    adata.obsm['X_scvi_denoised'] = full.sequential().get_sample_scale()
    adata.obsm['X_scvi_sample_rate'] = full.sequential().imputation()

    if return_posterior:
        return adata, full
    else:
        return adata

