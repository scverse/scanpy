from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy import linalg as la
from scverse_misc import Deprecation, deprecated_arg

from .. import logging as logg
from .._compat import CSBase
from .._docs import doc_out, doc_ref_compat, doc_use
from .._settings import Default
from .._utils import _doc_params, dim_acc, sanitize_anndata
from ..get._out import _legacy_out, _write_out
from ..get._use import _get_arr, _obs_df, _ref_arg, _resolve_obs, _resolve_ref

if TYPE_CHECKING:
    from collections.abc import Collection

    from anndata import AnnData

    from ..get._out import Prefix
    from ..get._types import AdRef, RepAcc


def _design_matrix(
    model: pd.DataFrame, batch_key: str, batch_levels: Collection[str]
) -> pd.DataFrame:
    """Compute a simple design matrix.

    Parameters
    ----------
    model
        Contains the batch annotation
    batch_key
        Name of the batch column
    batch_levels
        Levels of the batch annotation

    Returns
    -------
    The design matrix for the regression problem

    """
    import patsy

    design = patsy.dmatrix(
        f"~ 0 + C(Q('{batch_key}'), levels=batch_levels)",
        model,
        return_type="dataframe",
    )
    model = model.drop([batch_key], axis=1)
    numerical_covariates = model.select_dtypes("number").columns.array

    logg.info(f"Found {design.shape[1]} batches\n")
    other_cols = [c for c in model.columns.array if c not in numerical_covariates]

    if other_cols:
        col_repr = " + ".join(f"Q('{x}')" for x in other_cols)
        factor_matrix = patsy.dmatrix(
            f"~ {col_repr}", model[other_cols], return_type="dataframe"
        ).drop(columns="Intercept")

        design = pd.concat((design, factor_matrix), axis=1)
        logg.info(f"Found {len(other_cols)} categorical variables:")
        logg.info(f"\t{', '.join(other_cols)}\n")

    if numerical_covariates is not None:
        logg.info(f"Found {len(numerical_covariates)} numerical variables:")
        logg.info(f"\t{', '.join(numerical_covariates)}\n")

        for n_c in numerical_covariates:
            design[n_c] = model[n_c]

    return design


def _standardize_data(
    model: pd.DataFrame, data: pd.DataFrame, batch_key: str
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    """Standardize the data per gene.

    The aim here is to make mean and variance be comparable across batches.

    Parameters
    ----------
    model
        Contains the batch annotation
    data
        Contains the Data
    batch_key
        Name of the batch column in the model matrix

    Returns
    -------
    s_data
        Standardized Data
    design
        Batch assignment as one-hot encodings
    var_pooled
        Pooled variance per gene
    stand_mean
        Gene-wise mean

    """
    # compute the design matrix
    batch_items = model.groupby(batch_key, observed=True).groups.items()
    batch_levels, batch_info = zip(*batch_items, strict=True)
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))

    design = _design_matrix(model, batch_key, batch_levels)

    # compute pooled variance estimator
    b_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T)
    grand_mean = np.dot((n_batches / n_array).T, b_hat[:n_batch, :])
    var_pooled = (
        (data - np.dot(design, b_hat).T).pow(2).to_numpy().mean(axis=1, keepdims=True)
    )

    # Compute the means
    if np.sum(var_pooled == 0) > 0:
        print(f"Found {np.sum(var_pooled == 0)} genes with zero variance.")
    tmp = design.to_numpy(copy=True)
    tmp[:, :n_batch] = 0
    stand_mean = grand_mean[:, np.newaxis] + np.dot(tmp, b_hat).T

    # need to be a bit careful with the zero variance genes
    # just set the zero variance genes to zero in the standardized data
    s_data = np.where(
        var_pooled == 0,
        0,
        (data - stand_mean) / np.sqrt(var_pooled),
    )
    s_data = pd.DataFrame(s_data, index=data.index, columns=data.columns)

    return s_data, design, var_pooled, stand_mean


@_doc_params(
    use=doc_use("Which matrix to correct.", legacy=()),
    out=doc_out("combat"),
    ref=doc_ref_compat,
)
@deprecated_arg(
    "inplace", Deprecation("1.13.0", "Use `out=True` instead of `inplace=False`.")
)
@deprecated_arg("key", Deprecation("1.13.0", "Use `batch` instead."))
def combat(  # noqa: PLR0915
    adata: AnnData,
    batch: AdRef | str | Default = Default("A.obs['batch']"),
    *,
    covariates: Collection[AdRef | str] | None = None,
    use: RepAcc | str | None = None,
    out: RepAcc | str | Prefix | bool = False,
    # deprecated
    inplace: bool = True,
    key: str | None = None,
) -> np.ndarray | None:
    """ComBat function for batch effect correction :cite:p:`Johnson2006,Leek2012,Pedersen2012`.

    Corrects for batch effects by fitting linear models, gains statistical power
    via an EB framework where information is borrowed across genes.
    This uses the implementation `combat.py`_ :cite:p:`Pedersen2012`.

    .. _combat.py: https://github.com/brentp/combat.py

    .. array-support:: pp.combat

    Parameters
    ----------
    adata
        Annotated data matrix
    batch
        A categorical annotation that will be used for batch effect removal,
        e.g. `A.obs['batch']`. {ref}
    covariates
        Additional covariates besides the batch variable such as adjustment
        variables or biological condition (references or names like `batch`). This parameter refers to the design
        matrix `X` in Equation 2.1 in :cite:t:`Johnson2006` and to the `mod` argument in
        the original combat function in the sva R package.
        Note that not including covariates may introduce bias or lead to the
        removal of biological signal in unbalanced designs.
    {use}
    {out}
    inplace
        Whether to replace adata.X or to return the corrected data
    key
        The name of a `.obs` column, i.e. `key='b'` is `batch=A.obs['b']`.

    Returns
    -------
    Returns the corrected matrix (a :class:`numpy.ndarray`) if `out=True`, else `None`.

    """
    out = _legacy_out(out, "inplace", None if inplace else True)
    use = _resolve_obs(use)
    # check the input
    if isinstance(batch, Default):
        batch = None if key is not None else dim_acc("batch", dim="obs")
    # `batch` was positional before being renamed from `key`, so keep its strings’ meaning
    batch = _ref_arg(batch, key, name="batch", legacy_name="key", compat=True)
    covariates = [] if covariates is None else list(covariates)
    batch_ref, *cov_refs = (_resolve_ref(r) for r in [batch, *covariates])
    if batch_ref in cov_refs:
        msg = "Batch key and covariates cannot overlap"
        raise ValueError(msg)
    if len(cov_refs) != len(set(cov_refs)):
        msg = "Covariates must be unique"
        raise ValueError(msg)

    # only works on dense matrices so far
    x = _get_arr(adata, use)
    x = x.toarray().T if isinstance(x, CSBase) else x.T
    data = pd.DataFrame(data=x, columns=adata.obs_names)

    sanitize_anndata(adata)

    # construct a pandas series of the batch annotation
    try:
        model = _obs_df(adata, [batch, *covariates])
    except KeyError as e:
        msg = f"Could not find the batch or covariate {e} in adata.obs"
        raise ValueError(msg) from e
    # the display labels may not be valid in patsy formulas, so use positional ones
    batch_name, *_ = model.columns = [f"_{i}" for i in range(model.shape[1])]
    batch_info = model.groupby(batch_name, observed=True).indices
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info.values()])

    # check for batches with fewer than 2 cells
    small_batches = [
        batch for batch, size in zip(batch_info, n_batches, strict=True) if size < 2
    ]
    if small_batches:
        msg = (
            f"Batches {small_batches!r} have fewer than 2 cells. "
            "ComBat requires at least 2 cells per batch to estimate "
            "within-batch variance. Filter these batches before running combat."
        )
        raise ValueError(msg)

    # standardize across genes using a pooled variance estimator
    logg.info("Standardizing Data across genes.\n")
    s_data, design, var_pooled, stand_mean = _standardize_data(model, data, batch_name)

    # fitting the parameters on the standardized data
    logg.info("Fitting L/S model and finding priors\n")
    batch_design = design[design.columns[:n_batch]]
    # first estimate of the additive batch effect
    gamma_hat = (
        la.inv(batch_design.T @ batch_design) @ batch_design.T @ s_data.T
    ).to_numpy()
    # first estimate for the multiplicative batch effect
    delta_hat = [
        s_data.iloc[:, batch_idxs].var(axis=1) for batch_idxs in batch_info.values()
    ]

    # empirically fix the prior hyperparameters
    gamma_bar = gamma_hat.mean(axis=1)
    t2 = gamma_hat.var(axis=1)
    # a_prior and b_prior are the priors on lambda and theta from Johnson and Li (2006)
    a_prior = list(map(_aprior, delta_hat))
    b_prior = list(map(_bprior, delta_hat))

    logg.info("Finding parametric adjustments\n")
    # gamma star and delta star will be our empirical bayes (EB) estimators
    # for the additive and multiplicative batch effect per batch and cell
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_info.values()):
        # temp stores our estimates for the batch effect parameters.
        # temp[0] is the additive batch effect
        # temp[1] is the multiplicative batch effect
        gamma, delta = _it_sol(
            s_data.iloc[:, batch_idxs].values,
            gamma_hat[i],
            delta_hat[i].values,
            g_bar=gamma_bar[i],
            t2=t2[i],
            a=a_prior[i],
            b=b_prior[i],
        )

        gamma_star.append(gamma)
        delta_star.append(delta)

    logg.info("Adjusting data\n")
    bayesdata = s_data
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)

    # we now apply the parametric adjustment to the standardized data from above
    # loop over all batches in the data
    for j, batch_idxs in enumerate(batch_info.values()):
        # we basically subtract the additive batch effect, rescale by the ratio
        # of multiplicative batch effect to pooled variance and add the overall gene
        # wise mean
        dsq = np.sqrt(delta_star[j, :])
        numer = np.array(
            bayesdata.iloc[:, batch_idxs]
            - np.dot(batch_design.iloc[batch_idxs], gamma_star).T
        )
        bayesdata.iloc[:, batch_idxs] = numer / dsq[:, np.newaxis]

    bayesdata = bayesdata * np.sqrt(var_pooled) + stand_mean

    # put back into the adata object or return
    x = bayesdata.to_numpy().transpose()
    return _write_out(adata, x, out=out, base="combat", use=use)


def _it_sol(
    s_data: np.ndarray,
    g_hat: np.ndarray,
    d_hat: np.ndarray,
    *,
    g_bar: float,
    t2: float,
    a: float,
    b: float,
    conv: float = 0.0001,
) -> tuple[np.ndarray, np.ndarray]:
    """Iteratively compute the conditional posterior means for gamma and delta.

    gamma is an estimator for the additive batch effect, deltat is an estimator
    for the multiplicative batch effect. We use an EB framework to estimate these
    two. Analytical expressions exist for both parameters, which however depend on each other.
    We therefore iteratively evalutate these two expressions until convergence is reached.

    Parameters
    ----------
    s_data
        Contains the standardized Data
    g_hat
        Initial guess for gamma
    d_hat
        Initial guess for delta
    g_bar, t2, a, b
        Hyperparameters
    conv: float, optional (default: `0.0001`)
        convergence criterium

    Returns
    -------
    gamma
        estimated value for gamma
    delta
        estimated value for delta

    """  # noqa: D401
    n = (1 - np.isnan(s_data)).sum(axis=1)
    g_old = g_hat.copy()
    d_old = d_hat.copy()

    change = 1
    count = 0

    # They need to be initialized for numba to properly infer types
    g_new = g_old
    d_new = d_old
    # we place a normally distributed prior on gamma and and inverse gamma prior on delta
    # in the loop, gamma and delta are updated together. they depend on each other. we iterate until convergence.
    while change > conv:
        g_new = (t2 * n * g_hat + d_old * g_bar) / (t2 * n + d_old)
        sum2 = ((s_data - g_new[:, np.newaxis]) ** 2).sum(axis=1)
        d_new = (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)

        change = max(
            (abs(g_new - g_old) / g_old).max(), (abs(d_new - d_old) / d_old).max()
        )
        g_old = g_new  # .copy()
        d_old = d_new  # .copy()
        count = count + 1

    return g_new, d_new


def _aprior(delta_hat):
    m = delta_hat.mean()
    s2 = delta_hat.var()
    return (2 * s2 + m**2) / s2


def _bprior(delta_hat):
    m = delta_hat.mean()
    s2 = delta_hat.var()
    return (m * s2 + m**3) / s2
