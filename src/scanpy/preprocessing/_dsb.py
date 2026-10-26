"""dsb normalization of CITE-seq protein (ADT) counts :cite:p:`Mule2022`."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING, Literal, NamedTuple

import numba
import numpy as np
from fast_array_utils.conv import to_dense
from fast_array_utils.numba import njit

from .._compat import CSBase, warn
from .._utils import get_literal_vals, view_to_actual
from ..get import _get_arr, _set_obs_rep

if TYPE_CHECKING:
    from collections.abc import Iterable

    from anndata import AnnData
    from numpy.typing import DTypeLike, NDArray

    from .._utils.random import RNGLike, SeedLike

type _ScaleFactor = Literal["standardize", "mean_subtract"]

_SUBSET = 2000  # longer vectors start from a random subset this size
# EM iteration cap; only a guard against fits that never converge
_MAX_ITER = 100_000
_CHUNK = 65_536  # rows of empty droplets densified at a time
_ISOTYPE_PATTERN = re.compile(r"sotype|Iso|iso|control|CTRL|ctrl|Ctrl|ontrol")
_VAR_KEYS = ("dsb_ambient_mean", "dsb_ambient_sd")
_OBS_KEYS = ("dsb_cell_background_mean", "dsb_technical_component")


def dsb(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    adata_empty: AnnData | None = None,
    *,
    isotype_controls: str | Iterable[str] | None = None,
    use_isotype_controls: bool = True,
    denoise_counts: bool = True,
    pseudocount: float | None = None,
    scale_factor: _ScaleFactor = "standardize",
    fast_km: bool = False,
    quantile_clip: tuple[float, float] | None = None,
    layer: str | None = None,
    key_added: str | None = None,
    dtype: DTypeLike = np.float64,
    rng: SeedLike | RNGLike | None = None,
    copy: bool = False,
) -> AnnData | None:
    """Normalize and denoise protein (ADT) counts with dsb :cite:p:`Mule2022`.

    With empty droplets (`adata_empty`) this mirrors ``dsb::DSBNormalizeProtein``:

    1. *Ambient correction*: ``log(counts + pseudocount)`` of every protein is
       standardized with the mean and standard deviation of that protein in
       the empty droplets (or only centered, `scale_factor='mean_subtract'`).
    2. *Technical denoising* (`denoise_counts`): the first-component mean of a
       two-component Gaussian mixture fit to each cell's values is its
       background level. With the isotype controls it defines a technical
       component (first principal component of the standardized isotype values
       and background means), which is regressed out of every protein.

    Without `adata_empty` this mirrors ``dsb::ModelNegativeADTnorm``: step 1
    subtracts, per protein, the first-component mean of a mixture fit across
    cells.

    .. array-support:: pp.dsb

    Parameters
    ----------
    adata
        Cells-by-proteins raw ADT counts of cell-containing droplets.
    adata_empty
        Raw ADT counts of empty droplets with the same proteins. If `None`,
        the background is estimated from the cells.
    isotype_controls
        Name(s) of isotype control proteins in `adata.var_names` for step 2.
    use_isotype_controls
        Use `isotype_controls` in step 2, as recommended by dsb. If `False`, only
        the per-cell background mean is regressed out.
    denoise_counts
        Run step 2.
    pseudocount
        Added before log transformation. Defaults to 10 with `adata_empty` and
        1 without, as in dsb.
    scale_factor
        Step 1 with empty droplets: `'standardize'` or `'mean_subtract'`.
    fast_km
        Use the lower center of a two-cluster k-means as per-cell background.
        The exact optimum of the one-dimensional problem is computed.
    quantile_clip
        Clip every protein to these lower and upper quantiles (type 7) of its
        normalized values; dsb's `quantile.clipping` uses `(0.001, 0.9995)`.
    layer
        Layer with the counts in `adata` and `adata_empty`; `.X` if `None`.
    key_added
        Store the result in `adata.layers[key_added]` instead of overwriting
        `layer` (or `.X`).
    dtype
        Floating point type of the output. Computations use float64 and match
        R to ~1e-10.

        Where R would fail on a model fit, scanpy warns and continues: a cell
        whose mixture fit is singular in both models (R errors) gets the lower
        center of a two-cluster k-means as background, and EM stops after
        100,000 iterations (R has no limit).
    rng
        Random number generator for the mixture fits of vectors longer than
        2000 values, which start from a random subset. Without `adata_empty`
        and with more than 2000 cells, results vary between calls unless set.
    copy
        Return a modified copy instead of updating `adata` in place.

    Returns
    -------
    Returns `adata` if `copy=True`, else updates it and sets the following.
    Statistics of earlier calls are removed.

    `.X` or `.layers[layer | key_added]`
        dsb normalized values (dense, `dtype`).
    `.var['dsb_ambient_mean']`
        Per-protein background mean of ``log(counts + pseudocount)``: empty
        droplet mean, or first mixture component without `adata_empty` (0 for
        proteins whose fit failed).
    `.var['dsb_ambient_sd']`
        Empty droplet standard deviation (only with `adata_empty`).
    `.obs['dsb_cell_background_mean']`
        Per-cell background mean of step 2 (only with `denoise_counts`).
    `.obs['dsb_technical_component']`
        Covariate regressed out in step 2 (only with `denoise_counts`). Its sign
        is arbitrary in R; here it increases with the background mean.
    """
    if scale_factor not in (vals := get_literal_vals(_ScaleFactor)):
        msg = f"`scale_factor` must be one of {list(vals)}, got {scale_factor!r}."
        raise ValueError(msg)
    if adata_empty is None and scale_factor != "standardize":
        msg = "`scale_factor` only applies with `adata_empty`."
        raise ValueError(msg)
    if fast_km and not denoise_counts:
        msg = "`fast_km` only applies with `denoise_counts=True`."
        raise ValueError(msg)
    dtype = np.dtype(dtype)
    if dtype not in (np.float32, np.float64):
        msg = "`dtype` must be float32 or float64."
        raise ValueError(msg)
    if pseudocount is None:
        pseudocount = 1 if adata_empty is None else 10
    if not pseudocount > 0:
        msg = f"`pseudocount` must be positive, got {pseudocount}."
        raise ValueError(msg)
    if quantile_clip is not None and not 0 <= quantile_clip[0] < quantile_clip[1] <= 1:
        msg = f"`quantile_clip` must be increasing quantiles in [0, 1], got {quantile_clip}."
        raise ValueError(msg)
    proteins = adata.var_names
    if not proteins.is_unique:
        msg = (
            f"`adata.var_names` are not unique: {list(proteins[proteins.duplicated()])}"
        )
        raise ValueError(msg)
    iso_idx = _isotype_indices(
        proteins,
        isotype_controls,
        use_isotype_controls=use_isotype_controls,
        required=denoise_counts,
    )
    if adata.n_obs < 2 or (denoise_counts and adata.n_vars < 2):
        msg = "dsb needs at least 2 cells, and 2 proteins with `denoise_counts`."
        raise ValueError(msg)
    if adata_empty is not None:
        adata_empty = _match_proteins(adata_empty, proteins)
        if adata_empty.n_obs < 2:
            msg = "`adata_empty` needs at least 2 empty droplets."
            raise ValueError(msg)
        if layer is not None and layer not in adata_empty.layers:
            msg = f"`adata_empty` has no layer {layer!r}."
            raise ValueError(msg)
    rng = np.random.default_rng(rng)

    if copy:
        adata = adata.copy()
    else:
        view_to_actual(adata)
    adata.var.drop(columns=[k for k in _VAR_KEYS if k in adata.var], inplace=True)
    adata.obs.drop(columns=[k for k in _OBS_KEYS if k in adata.obs], inplace=True)
    norm = _log_counts(_get_arr(adata, layer=layer), pseudocount, "adata")

    if adata_empty is not None:
        mean, sd, constant = _empty_stats(
            _get_arr(adata_empty, layer=layer), pseudocount
        )
        var_stats = {"dsb_ambient_mean": mean, "dsb_ambient_sd": sd}
        norm -= mean
        if scale_factor == "standardize":
            if constant.any():
                msg = (
                    "Proteins with the same count in every empty droplet cannot be "
                    f"standardized; remove them from both objects: {list(proteins[constant])}"
                )
                raise ValueError(msg)
            if (low := sd < 0.05).any():
                warn(
                    "Proteins with low background log-count sd < 0.05, check "
                    f"their raw and normalized distributions: {list(proteins[low])}",
                    UserWarning,
                )
            norm /= sd
    else:
        mean = _mixture_fits(norm.T, rng).params[:, 0]
        if (failed := np.isnan(mean)).any():
            warn(
                f"The background could not be fit for proteins "
                f"{list(proteins[failed])}; they are only log transformed.",
                UserWarning,
            )
            mean[failed] = 0
        var_stats = {"dsb_ambient_mean": mean}
        norm -= mean

    obs_stats = {}
    if denoise_counts:
        if fast_km:
            background = _kmeans_low_centers(norm)
        else:
            background = _mixture_fits(norm, rng).params[:, 0]
            if (failed := np.isnan(background)).any():
                warn(
                    f"The mixture model could not be fit for {int(failed.sum())} "
                    "cells; using the lower k-means center of their values instead.",
                    UserWarning,
                )
                background[failed] = _kmeans_low_centers(norm[failed])
        technical = background
        if iso_idx is not None:
            noise = np.concatenate([norm[:, iso_idx], background[:, None]], axis=1)
            if (constant := np.ptp(noise, axis=0) == 0).any():
                names = [*proteins[iso_idx], "cell background mean"]
                msg = (
                    "Cannot standardize the technical component, constant across "
                    f"cells: {[n for n, c in zip(names, constant, strict=True) if c]}"
                )
                raise ValueError(msg)
            technical = _first_pc_scores(noise)
        # subtract the slope times the centered technical component;
        # a constant covariate leaves the data unchanged
        if np.ptp(technical) > 0:
            vc = technical - technical.mean()
            slopes = (vc @ norm) / (vc @ vc)
            for start in range(0, norm.shape[0], _CHUNK):
                rows = slice(start, start + _CHUNK)
                norm[rows] -= vc[rows, None] * slopes
        else:
            warn(
                "The technical component is constant; step 2 leaves the data "
                "unchanged.",
                UserWarning,
            )
        obs_stats = {
            "dsb_cell_background_mean": background,
            "dsb_technical_component": technical,
        }

    if quantile_clip is not None:
        q = np.quantile(norm, quantile_clip, axis=0)
        np.clip(norm, q[0], q[1], out=norm)

    out = norm.astype(dtype, copy=False)
    if key_added is None:
        _set_obs_rep(adata, out, layer=layer)
    else:
        adata.layers[key_added] = out
    for key, val in var_stats.items():
        adata.var[key] = val.astype(dtype, copy=False)
    for key, val in obs_stats.items():
        adata.obs[key] = val.astype(dtype, copy=False)
    return adata if copy else None


def _isotype_indices(
    proteins, isotype_controls, *, use_isotype_controls: bool, required: bool
) -> NDArray[np.intp] | None:
    if isinstance(isotype_controls, str):
        isotype_controls = [isotype_controls]
    elif isotype_controls is not None:
        isotype_controls = list(isotype_controls)
    if isotype_controls:
        if not use_isotype_controls:
            msg = "Pass `isotype_controls` or `use_isotype_controls=False`, not both."
            raise ValueError(msg)
        if missing := [p for p in isotype_controls if p not in proteins]:
            msg = f"isotype controls not found in `adata.var_names`: {missing}"
            raise ValueError(msg)
        return proteins.get_indexer(isotype_controls)
    if required and use_isotype_controls:
        detected = [p for p in proteins if _ISOTYPE_PATTERN.search(p)]
        msg = (
            "Set `isotype_controls` to the names of the isotype controls, or "
            "`use_isotype_controls=False` to regress out only the per-cell "
            f"background (not recommended). Potential isotype controls: {detected}"
        )
        raise ValueError(msg)
    return None


def _log_counts(x, pseudocount: float, name: str) -> NDArray[np.float64]:
    """``log(x + pseudocount)`` as a new dense float64 array."""
    if x is None:
        msg = f"`{name}.X` is None; pass the `layer` with the counts."
        raise ValueError(msg)
    # only densified sparse matrices are owned; never modify the input in place
    owned = isinstance(x, CSBase)
    out = np.array(
        to_dense(x, order="C", to_cpu_memory=True),
        dtype=np.float64,
        order="C",
        copy=None if owned else True,
    )
    if not (np.isfinite(out).all() and (out >= 0).all()):
        msg = f"`{name}` must contain raw counts: finite and non-negative values."
        raise ValueError(msg)
    out += pseudocount
    return np.log(out, out=out)


def _empty_stats(
    x, pseudocount: float
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.bool]]:
    """Per-protein mean, sd and constancy of the log counts, in row chunks."""
    n, mean, m2 = 0, 0.0, 0.0
    lo, hi = np.inf, -np.inf
    for start in range(0, x.shape[0], _CHUNK):
        chunk = _log_counts(x[start : start + _CHUNK], pseudocount, "adata_empty")
        n_b, mean_b = chunk.shape[0], chunk.mean(axis=0)
        m2_b = ((chunk - mean_b) ** 2).sum(axis=0)
        # merge with the previous chunks (Chan et al.)
        delta, n_ab = mean_b - mean, n + n_b
        mean = mean + delta * (n_b / n_ab)
        m2 = m2 + m2_b + delta**2 * (n * n_b / n_ab)
        n = n_ab
        lo, hi = np.minimum(lo, chunk.min(axis=0)), np.maximum(hi, chunk.max(axis=0))
    return mean, np.sqrt(m2 / (n - 1)), lo == hi


def _match_proteins(adata_empty: AnnData, proteins) -> AnnData:
    names = adata_empty.var_names
    if names.equals(proteins):
        return adata_empty
    if (
        not names.is_unique
        or len(names) != len(proteins)
        or set(names) != set(proteins)
    ):
        diff = sorted(set(names).symmetric_difference(proteins))
        dup = list(names[names.duplicated()])
        msg = (
            "`adata` and `adata_empty` have different proteins: "
            f"{diff[:20]}{f', duplicated: {dup[:20]}' if dup else ''}"
        )
        raise ValueError(msg)
    warn("Reordered the proteins of `adata_empty` to match `adata`.", UserWarning)
    return adata_empty[:, proteins]


class _MixtureFit(NamedTuple):
    params: NDArray[
        np.float64
    ]  # (n, 6): mean0, mean1, var0, var1, pro0, pro1; NaN = no fit
    capped: NDArray[np.bool]  # EM of either model hit the iteration cap


def _mixture_fits(
    m: NDArray[np.float64], rng: np.random.Generator, *, max_iter: int = _MAX_ITER
) -> _MixtureFit:
    """Fit a univariate two-component Gaussian mixture to every row of `m`.

    * Start: hard split ``x < q``, ``q`` the type 7 quantile at the fraction
      ``j / k`` with the smallest ``k >= 2`` for which ``min < q < max`` (the
      median for untied data). ``mean[1]`` is the component started on the
      lower class.
    * EM (M-step, E-step) until ``|l_k - l_{k-1}| / (1 + |l_k|) < 1e-5``; a
      model with a variance ``<= eps`` is dropped. One more M-step is applied
      when ``sum((pro - colMeans(z))^2) > sqrt(eps)``. The larger
      ``2 loglik - npar log(n)`` of models E and V wins.
    * Vectors longer than 2000 values: split and first M-step on a random
      subset of 2000, EM on all values, extra M-step always applied.

    EM runs once per distinct value, weighted by its multiplicity.
    """
    n_fits, n = m.shape
    if n > _SUBSET:
        sub = np.stack([rng.choice(n, _SUBSET, replace=False) for _ in range(n_fits)])
    else:
        sub = np.empty((n_fits, 0), dtype=np.intp)
    params, capped = _mix2_fits(m, sub, max_iter)
    if n_capped := int(capped.sum()):
        warn(
            f"EM stopped at the iteration limit ({max_iter}) without "
            f"converging for {n_capped} fits; these may differ from R.",
            UserWarning,
        )
    return _MixtureFit(params, capped)


# Mixture parameters are passed around as float64 tuples
# (mu0, mu1, var0, var1, pro0, pro1).
_TOL = 1e-5  # relative log-likelihood change
_EPS = float(np.finfo(np.float64).eps)  # variances at or below are singular
# One more M-step when sum((pro - colMeans(z))^2) exceeds sqrt(DBL_EPSILON).
_EXTRA_MSTEP_TOL = math.sqrt(_EPS)


@njit
def _mix2_fits(
    m: NDArray[np.float64], sub: NDArray[np.intp], max_iter: int
) -> tuple[NDArray[np.float64], NDArray[np.bool]]:
    n_fits, n = m.shape
    subset = sub.shape[1] > 0
    params = np.full((n_fits, 6), np.nan)
    capped = np.zeros(n_fits, dtype=np.bool)
    for f in numba.prange(n_fits):
        xs = np.sort(m[f])
        values, weights = _unique_counts(xs)
        start = np.sort(m[f][sub[f]]) if subset else xs
        q = _threshold(start)
        best_bic = -np.inf
        for model in range(2):
            p, ll, it = _em(
                values,
                weights,
                _start_mstep(start, q, model == 0),
                model == 0,
                subset,
                max_iter,
            )
            capped[f] |= it >= max_iter
            bic = 2 * ll - (4 + model) * math.log(n)
            if it >= 0 and bic > best_bic:
                best_bic = bic
                for k in range(6):
                    params[f, k] = p[k]
    return params, capped


@numba.njit(cache=True)  # noqa: TID251
def _unique_counts(xs: NDArray[np.float64]):
    """Distinct values of the sorted `xs` and their multiplicities."""
    values = np.empty_like(xs)
    weights = np.empty_like(xs)
    k = -1
    for i in range(xs.shape[0]):
        if k < 0 or xs[i] != values[k]:
            k += 1
            values[k] = xs[i]
            weights[k] = 0.0
        weights[k] += 1.0
    return values[: k + 1], weights[: k + 1]


@numba.njit(cache=True)  # noqa: TID251
def _threshold(xs: NDArray[np.float64]) -> float:
    """Start split of the sorted values `xs`.

    The type 7 quantile at the fraction ``j / k`` with the smallest ``k >= 2``
    for which ``min < q < max`` (the median for untied data); NaN when all
    values are equal.
    """
    n = xs.shape[0]
    if n < 2 or not xs[0] < xs[n - 1]:
        return np.nan
    a = 1  # copies of the min
    while xs[a] == xs[0]:
        a += 1
    b = n - 1  # first index of the max
    while xs[b - 1] == xs[n - 1]:
        b -= 1
    # smallest k with an integer j in ((a - 1) k / (n - 1), b k / (n - 1))
    k = 2
    j = (a - 1) * k // (n - 1) + 1
    while j * (n - 1) >= b * k:
        k += 1
        j = (a - 1) * k // (n - 1) + 1
    lo, rem = (n - 1) * j // k, (n - 1) * j % k
    if rem == 0:
        return xs[lo]
    return xs[lo] + rem / k * (xs[lo + 1] - xs[lo])


@numba.njit(cache=True, error_model="numpy")  # noqa: TID251
def _finish_mstep(w0, w1, a0, a1, b0, b1, n, equal, sh0, sh1):  # noqa: PLR0917
    """M-step from moments about shifts sh_k.

    ``w_k = sum z_k``, ``a_k = sum z_k (x - sh_k)``, ``b_k = sum z_k (x - sh_k)^2``.
    """
    d0, d1 = a0 / w0, a1 / w1
    ss0, ss1 = b0 - a0 * d0, b1 - a1 * d1
    s0 = (ss0 + ss1) / n if equal else ss0 / w0
    s1 = (ss0 + ss1) / n if equal else ss1 / w1
    return (sh0 + d0, sh1 + d1, s0, s1, w0 / n, w1 / n)


@numba.njit(cache=True, error_model="numpy")  # noqa: TID251
def _start_mstep(x: NDArray[np.float64], q: float, equal: bool):  # noqa: FBT001
    """M-step from the hard split ``x < q``."""
    w0 = w1 = a0 = a1 = b0 = b1 = 0.0
    for i in range(x.shape[0]):
        d = x[i] - q
        if x[i] < q:
            w0 += 1.0
            a0 += d
            b0 += d * d
        else:
            w1 += 1.0
            a1 += d
            b1 += d * d
    return _finish_mstep(w0, w1, a0, a1, b0, b1, x.shape[0], equal, q, q)


@numba.njit(cache=True)  # noqa: TID251
def _valid(p) -> bool:
    finite = (
        math.isfinite(p[0])
        and math.isfinite(p[1])
        and math.isfinite(p[2])
        and math.isfinite(p[3])
    )
    return finite and p[2] > _EPS and p[3] > _EPS and p[4] > 0 and p[5] > 0


@numba.njit(cache=True, error_model="numpy")  # noqa: TID251
def _em(values, weights, p, equal, subset, max_iter):  # noqa: PLR0917
    """Run EM from `p` on `values` with multiplicities `weights`.

    Returns the parameters, log-likelihood and iterations (-1 if singular or
    degenerate). With `subset`, `p` is the start of a subset fit, whose
    log-likelihood is not compared and whose extra M-step is always applied.
    """
    if not _valid(p):
        return p, np.nan, -1
    n = weights.sum()
    ll_prev, it = 0.0, -1 if subset else 0
    while True:
        # E-step at p fused with the next M-step: moments about the current means
        c0 = math.log(p[4]) - 0.5 * math.log(2 * math.pi * p[2])
        c1 = math.log(p[5]) - 0.5 * math.log(2 * math.pi * p[3])
        h0, h1 = 0.5 / p[2], 0.5 / p[3]
        w0 = w1 = a0 = a1 = b0 = b1 = ll = 0.0
        for i in range(values.shape[0]):
            d0, d1 = values[i] - p[0], values[i] - p[1]
            l0, l1 = c0 - d0 * d0 * h0, c1 - d1 * d1 * h1
            first = l0 > l1
            hi, lo = (l0, l1) if first else (l1, l0)
            e = math.exp(lo - hi)
            inv = weights[i] / (1.0 + e)
            z0, z1 = (inv, e * inv) if first else (e * inv, inv)
            ll += weights[i] * (hi + math.log1p(e))
            w0 += z0
            w1 += z1
            a0 += z0 * d0
            a1 += z1 * d1
            b0 += z0 * d0 * d0
            b1 += z1 * d1 * d1
        pn = _finish_mstep(w0, w1, a0, a1, b0, b1, n, equal, p[0], p[1])
        if not math.isfinite(ll):
            return p, np.nan, -1
        it += 1
        if it > 0:
            if (it > 1 and abs(ll - ll_prev) < _TOL * (1 + abs(ll))) or it >= max_iter:
                # pn's proportions are the column means of the final responsibilities
                dp = (p[4] - pn[4]) ** 2 + (p[5] - pn[5]) ** 2
                if (subset or dp > _EXTRA_MSTEP_TOL) and _valid(pn):
                    p = pn
                return p, ll, it
            ll_prev = ll
        if not _valid(pn):
            return p, np.nan, -1
        p = pn


@njit
def _kmeans_low_centers(m: NDArray[np.float64]) -> NDArray[np.float64]:
    """Lower center of the optimal 1-D two-means split of every row."""
    n_rows, n = m.shape
    out = np.empty(n_rows)
    for r in numba.prange(n_rows):
        s = np.sort(m[r])
        total = s.sum()
        left, best, low = 0.0, -np.inf, s[0]
        for k in range(1, n):
            left += s[k - 1]
            # minimizing the within-cluster sum of squares maximizes this;
            # splits inside runs of tied values are not real partitions
            if s[k - 1] < s[k]:
                score = left * left / k + (total - left) ** 2 / (n - k)
                if score > best:
                    best, low = score, left / k
        out[r] = low
    return out


def _first_pc_scores(a: NDArray[np.float64]) -> NDArray[np.float64]:
    """First principal component scores of the standardized columns.

    Oriented along the last column (the background mean).
    """
    z = (a - a.mean(axis=0)) / a.std(axis=0, ddof=1)
    v = np.linalg.eigh(z.T @ z)[1][:, -1]
    return z @ (v if v[-1] >= 0 else -v)
