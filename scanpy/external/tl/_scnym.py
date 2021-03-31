"""
Classify cell identities using scNym
"""


def scnym(*args, **kwargs) -> None:
    """
    scNym: Semi-supervised adversarial neural networks for
    single cell classification [Kimmel2020]_.

    scNym is a cell identity classifier that transfers annotations from one
    single cell experiment to another. The model is implemented as a neural
    network that employs MixMatch semi-supervision and a domain adversary to
    take advantage of unlabeled data during training. scNym offers superior
    performance to many baseline single cell identity classification methods.

    Parameters
    ----------
    adata
        Annotated data matrix used for training or prediction.
        If `"scNym_split"` in `.obs_keys()`, uses the cells annotated
        `"train", "val"` to select data splits.
    task
        Task to perform, either "train" or "predict".
        If "train", uses `adata` as labeled training data.
        If "predict", uses `trained_model` to infer cell identities for
        observations in `adata`.
    groupby
        Column in `adata.obs` that contains cell identity annotations.
        Values of `"Unlabeled"` indicate that a given cell should be used
        only as unlabeled data during training.
    domain_groupby
        Column in `adata.obs` that contains domain labels as integers.
        Each domain of origin (e.g. batch, species) should be given a unique
        domain label.
        If `domain_groupby is None`, train and target data are each considered
        a unique domain.
    out_path
        Path to a directory for saving scNym model weights and training logs.
    trained_model
        Path to the output directory of an scNym training run
        or a string specifying a pretrained model.
        If provided while `task == "train"`, used as an initialization.
    config
        Configuration name or dictionary of configuration of parameters.
        Pre-defined configurations:
            "new_identity_discovery" - Default. Employs pseudolabel thresholding to
            allow for discovery of new cell identities in the target dataset using
            scNym confidence scores.
            "no_new_identity" - Assumes all cells in the target data belong to one
            of the classes in the training data. Recommended to improve performance
            when this assumption is valid.
    key_added
        Key added to `adata.obs` with scNym predictions if `task=="predict"`.
    copy
        copy the AnnData object before predicting cell types.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    `X_scnym` : :class:`~numpy.ndarray`, (:attr:`~anndata.AnnData.obsm`, shape=(n_samples, n_hidden), dtype `float`)
        scNym embedding coordinates of data.
    `scNym` : (`adata.obs`, dtype `str`)
        scNym cell identity predictions for each observation.
    `scNym_train_results` : :class:`~dict`, (:attr:`~anndata.AnnData.uns`)
        results of scNym model training.

    Examples
    --------
    >>> import scanpy as sc
    >>> from scnym.api import scnym_api, atlas2target

    **Loading Data and preparing labels**

    >>> adata = sc.datasets.kang17()
    >>> target_bidx = adata.obs['stim']=='stim'
    >>> adata.obs['cell'] = np.array(adata.obs['cell'])
    >>> adata.obs.loc[target_bidx, 'cell'] = 'Unlabeled'

    **Train an scNym model**

    >>> scnym_api(
    ...   adata=adata,
    ...   task='train',
    ...   groupby='clusters',
    ...   out_path='./scnym_outputs',
    ...   config='no_new_identity',
    ... )

    **Predict cell identities with the trained scNym model**

    >>> path_to_model = './scnym_outputs/'
    >>> scnym_api(
    ...   adata=adata,
    ...   task='predict',
    ...   groupby='scNym',
    ...   trained_model=path_to_model,
    ...   config='no_new_identity',
    ... )

    **Perform semi-supervised training with an atlas**

    >>> joint_adata = atlas2target(
    ...   adata=adata,
    ...   species='mouse',
    ...   key_added='annotations',
    ... )
    >>> scnym_api(
    ...   adata=joint_adata,
    ...   task='train',
    ...   groupby='annotations',
    ...   out_path='./scnym_outputs',
    ...   config='no_new_identity',
    ... )
    """
    try:
        import scnym
    except ImportError:
        print("Please install `scnym`:\n\t`pip install scnym`.")

    # native `scnym` API is `scanpy` compatibile
    # modifies `adata` in place, no returns
    scnym.api.scnym_api(*args, **kwargs)
    return
