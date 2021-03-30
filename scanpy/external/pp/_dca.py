from types import MappingProxyType
from typing import Optional, Sequence, Union, Mapping, Any

from anndata import AnnData

from ..._compat import Literal
from ..._utils import AnyRandom


_AEType = Literal['zinb-conddisp', 'zinb', 'nb-conddisp', 'nb']


def dca(
    adata: AnnData,
    mode: Literal['denoise', 'latent'] = 'denoise',
    ae_type: _AEType = 'nb-conddisp',
    normalize_per_cell: bool = True,
    scale: bool = True,
    log1p: bool = True,
    # network args
    hidden_size: Sequence[int] = (64, 32, 64),
    hidden_dropout: Union[float, Sequence[float]] = 0.0,
    batchnorm: bool = True,
    activation: str = 'relu',
    init: str = 'glorot_uniform',
    network_kwds: Mapping[str, Any] = MappingProxyType({}),
    # training args
    epochs: int = 300,
    reduce_lr: int = 10,
    early_stop: int = 15,
    batch_size: int = 32,
    optimizer: str = 'RMSprop',
    random_state: AnyRandom = 0,
    threads: Optional[int] = None,
    learning_rate: Optional[float] = None,
    verbose: bool = False,
    training_kwds: Mapping[str, Any] = MappingProxyType({}),
    return_model: bool = False,
    return_info: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Deep count autoencoder [Eraslan18]_.

    Fits a count autoencoder to the raw count data given in the anndata object
    in order to denoise the data and to capture hidden representation of
    cells in low dimensions. Type of the autoencoder and return values are
    determined by the parameters.

    .. note::
        More information and bug reports `here <https://github.com/theislab/dca>`__.

    Parameters
    ----------
    adata
        An anndata file with `.raw` attribute representing raw counts.
    mode
        `denoise` overwrites `adata.X` with denoised expression values.
        In `latent` mode DCA adds `adata.obsm['X_dca']` to given adata
        object. This matrix represent latent representation of cells via DCA.
    ae_type
        Type of the autoencoder. Return values and the architecture is
        determined by the type e.g. `nb` does not provide dropout
        probabilities. Types that end with "-conddisp", assumes that dispersion is mean dependant.
    normalize_per_cell
        If true, library size normalization is performed using
        the `sc.pp.normalize_per_cell` function in Scanpy and saved into adata
        object. Mean layer is re-introduces library size differences by
        scaling the mean value of each cell in the output layer. See the
        manuscript for more details.
    scale
        If true, the input of the autoencoder is centered using
        `sc.pp.scale` function of Scanpy. Note that the output is kept as raw
        counts as loss functions are designed for the count data.
    log1p
        If true, the input of the autoencoder is log transformed with a
        pseudocount of one using `sc.pp.log1p` function of Scanpy.
    hidden_size
        Width of hidden layers.
    hidden_dropout
        Probability of weight dropout in the autoencoder (per layer if list
        or tuple).
    batchnorm
        If true, batch normalization is performed.
    activation
        Activation function of hidden layers.
    init
        Initialization method used to initialize weights.
    network_kwds
        Additional keyword arguments for the autoencoder.
    epochs
        Number of total epochs in training.
    reduce_lr
        Reduces learning rate if validation loss does not improve in given number of epochs.
    early_stop
        Stops training if validation loss does not improve in given number of epochs.
    batch_size
        Number of samples in the batch used for SGD.
    optimizer
        Type of optimization method used for training.
    random_state
        Seed for python, numpy and tensorflow.
    threads
        Number of threads to use in training. All cores are used by default.
    learning_rate
        Learning rate to use in the training.
    verbose
        If true, prints additional information about training and architecture.
    training_kwds
        Additional keyword arguments for the training process.
    return_model
        If true, trained autoencoder object is returned. See "Returns".
    return_info
        If true, all additional parameters of DCA are stored in `adata.obsm` such as dropout
        probabilities (obsm['X_dca_dropout']) and estimated dispersion values
        (obsm['X_dca_dispersion']), in case that autoencoder is of type
        zinb or zinb-conddisp.
    copy
        If true, a copy of anndata is returned.

    Returns
    -------
    If `copy` is true and `return_model` is false, AnnData object is returned.

    In "denoise" mode, `adata.X` is overwritten with the denoised values.
    In "latent" mode, latent low dimensional representation of cells are stored
    in `adata.obsm['X_dca']` and `adata.X` is not modified.
    Note that these values are not corrected for library size effects.

    If `return_info` is true, all estimated distribution parameters are stored
    in AnnData like this:

    `.obsm["X_dca_dropout"]`
        The mixture coefficient (pi) of the zero component in ZINB,
        i.e. dropout probability (if `ae_type` is `zinb` or `zinb-conddisp`).
    `.obsm["X_dca_dispersion"]`
        The dispersion parameter of NB.
    `.uns["dca_loss_history"]`
        The loss history of the training.
        See `.history` attribute of Keras History class for mode details.

    Finally, the raw counts are stored in `.raw` attribute of AnnData object.

    If `return_model` is given, trained model is returned.
    When both `copy` and `return_model` are true,
    a tuple of anndata and model is returned in that order.
    """

    try:
        from dca.api import dca
    except ImportError:
        raise ImportError('Please install dca package (>= 0.2.1) via `pip install dca`')

    return dca(
        adata,
        mode=mode,
        ae_type=ae_type,
        normalize_per_cell=normalize_per_cell,
        scale=scale,
        log1p=log1p,
        hidden_size=hidden_size,
        hidden_dropout=hidden_dropout,
        batchnorm=batchnorm,
        activation=activation,
        init=init,
        network_kwds=network_kwds,
        epochs=epochs,
        reduce_lr=reduce_lr,
        early_stop=early_stop,
        batch_size=batch_size,
        optimizer=optimizer,
        random_state=random_state,
        threads=threads,
        learning_rate=learning_rate,
        verbose=verbose,
        training_kwds=training_kwds,
        return_model=return_model,
        return_info=return_info,
        copy=copy,
    )
