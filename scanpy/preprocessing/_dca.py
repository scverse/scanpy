def dca(adata,
        mode='denoise',
        ae_type='zinb-conddisp',
        normalize_per_cell=True,
        scale=True,
        log1p=True,
        # network args
        hidden_size=(64, 32, 64),
        hidden_dropout=0.,
        batchnorm=True,
        activation='relu',
        init='glorot_uniform',
        network_kwds={},
        # training args
        epochs=300,
        reduce_lr=10,
        early_stop=15,
        batch_size=32,
        optimizer='rmsprop',
        random_state=0,
        threads=None,
        verbose=False,
        training_kwds={},
        return_model=False,
        return_info=False,
        copy=False
        ):
    """Deep count autoencoder [Eraslan18]_.

    Fits a count autoencoder to the raw count data given in the anndata object
    in order to denoise the data and to capture hidden representation of
    cells in low dimensions. Type of the autoencoder and return values are
    determined by the parameters.
    
    .. note::
        More information and bug reports `here <https://github.com/theislab/dca>`__.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        An anndata file with `.raw` attribute representing raw counts.
    mode : `str`, optional. `denoise`(default), or `latent`.
        `denoise` overwrites `adata.X` with denoised expression values.
        In `latent` mode DCA adds `adata.obsm['X_dca']` to given adata
        object. This matrix represent latent representation of cells via DCA.
    ae_type : `str`, optional. `zinb-conddisp`(default), `zinb`, `nb-conddisp` or `nb`.
        Type of the autoencoder. Return values and the architecture is
        determined by the type e.g. `nb` does not provide dropout
        probabilities. Types that end with "-conddisp", assumes that dispersion is mean dependant.
    normalize_per_cell : `bool`, optional. Default: `True`.
        If true, library size normalization is performed using
        the `sc.pp.normalize_per_cell` function in Scanpy and saved into adata
        object. Mean layer is re-introduces library size differences by
        scaling the mean value of each cell in the output layer. See the
        manuscript for more details.
    scale : `bool`, optional. Default: `True`.
        If true, the input of the autoencoder is centered using
        `sc.pp.scale` function of Scanpy. Note that the output is kept as raw
        counts as loss functions are designed for the count data.
    log1p : `bool`, optional. Default: `True`.
        If true, the input of the autoencoder is log transformed with a
        pseudocount of one using `sc.pp.log1p` function of Scanpy.
    hidden_size : `tuple` or `list`, optional. Default: (64, 32, 64).
        Width of hidden layers.
    hidden_dropout : `float`, `tuple` or `list`, optional. Default: 0.0.
        Probability of weight dropout in the autoencoder (per layer if list
        or tuple).
    batchnorm : `bool`, optional. Default: `True`.
        If true, batch normalization is performed.
    activation : `str`, optional. Default: `relu`.
        Activation function of hidden layers.
    init : `str`, optional. Default: `glorot_uniform`.
        Initialization method used to initialize weights.
    network_kwds : `dict`, optional.
        Additional keyword arguments for the autoencoder.
    epochs : `int`, optional. Default: 300.
        Number of total epochs in training.
    reduce_lr : `int`, optional. Default: 10.
        Reduces learning rate if validation loss does not improve in given number of epochs.
    early_stop : `int`, optional. Default: 15.
        Stops training if validation loss does not improve in given number of epochs.
    batch_size : `int`, optional. Default: 32.
        Number of samples in the batch used for SGD.
    optimizer : `str`, optional. Default: "rmsprop".
        Type of optimization method used for training.
    random_state : `int`, optional. Default: 0.
        Seed for python, numpy and tensorflow.
    threads : `int` or None, optional. Default: None
        Number of threads to use in training. All cores are used by default.
    verbose : `bool`, optional. Default: `False`.
        If true, prints additional information about training and architecture.
    training_kwds : `dict`, optional.
        Additional keyword arguments for the training process.
    return_model : `bool`, optional. Default: `False`.
        If true, trained autoencoder object is returned. See "Returns".
    return_info : `bool`, optional. Default: `False`.
        If true, all additional parameters of DCA are stored in `adata.obsm` such as dropout
        probabilities (obsm['X_dca_dropout']) and estimated dispersion values
        (obsm['X_dca_dispersion']), in case that autoencoder is of type
        zinb or zinb-conddisp.
    copy : `bool`, optional. Default: `False`.
        If true, a copy of anndata is returned.

    Returns
    -------
    If `copy` is true and `return_model` is false, AnnData object is returned.

    In "denoise" mode, `adata.X` is overwritten with the denoised values. In "latent" mode, latent\
    low dimensional representation of cells are stored in `adata.obsm['X_dca']` and `adata.X`\
    is not modified. Note that these values are not corrected for library size effects.

    If `return_info` is true, all estimated distribution parameters are stored in AnnData such as:

    - `.obsm["X_dca_dropout"]` which is the mixture coefficient (pi) of the zero component\
    in ZINB, i.e. dropout probability (only if `ae_type` is `zinb` or `zinb-conddisp`).

    - `.obsm["X_dca_dispersion"]` which is the dispersion parameter of NB.

    - `.uns["dca_loss_history"]` which stores the loss history of the training. See `.history`\
    attribute of Keras History class for mode details.

    Finally, the raw counts are stored in `.raw` attribute of AnnData object.

    If `return_model` is given, trained model is returned. When both `copy` and `return_model`\
    are true, a tuple of anndata and model is returned in that order.
    """

    try:
        from dca.api import dca
    except ImportError:
        raise ImportError('Please install dca package (>= 0.2.1) via `pip install dca`')

    return dca(adata,
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
        verbose=verbose,
        training_kwds=training_kwds,
        return_model=return_model)
