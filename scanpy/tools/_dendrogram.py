"""
Computes a dendrogram based on a given categorical observation.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal

import pandas as pd
from pandas.api.types import CategoricalDtype

from .. import logging as logg
from .._compat import old_positionals
from .._utils import _doc_params
from ..neighbors._doc import doc_n_pcs, doc_use_rep
from ._utils import _choose_representation, _resolve_axis

if TYPE_CHECKING:
    from anndata import AnnData


@old_positionals(
    "n_pcs",
    "use_rep",
    "var_names",
    "use_raw",
    "cor_method",
    "linkage_method",
    "optimal_ordering",
    "key_added",
    "inplace",
)
@_doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def dendrogram(
    adata: AnnData,
    groupby: str | Sequence[str] | None = None,
    *,
    axis: Literal["obs", 0, "var", 1] = "obs",
    n_pcs: int | None = None,
    use_rep: str | None = None,
    var_names: Sequence[str] | None = None,
    use_raw: bool | None = None,
    cor_method: str = "pearson",
    linkage_method: str = "complete",
    optimal_ordering: bool = False,
    key_added: str | None = None,
    inplace: bool = True,
) -> dict[str, Any] | None:
    """\
    Computes a hierarchical clustering for the given `groupby` categories.

    By default, the PCA representation is used unless `.X`
    has less than 50 variables.

    Alternatively, a list of `var_names` (e.g. genes) can be given.

    Average values of either `var_names` or components are used
    to compute a correlation matrix.

    The hierarchical clustering can be visualized using
    :func:`scanpy.pl.dendrogram` or multiple other visualizations that can
    include a dendrogram: :func:`~scanpy.pl.matrixplot`,
    :func:`~scanpy.pl.heatmap`, :func:`~scanpy.pl.dotplot`,
    and :func:`~scanpy.pl.stacked_violin`.

    .. note::
        The computation of the hierarchical clustering is based on predefined
        groups and not per cell. The correlation matrix is computed using by
        default pearson but other methods are available.

    Parameters
    ----------
    adata
        Annotated data matrix
    groupby
        The obs column(s) to use to group observations. Default is None.
    axis
        Axis along which to calculate the dendrogram.
    {n_pcs}
    {use_rep}
    var_names
        List of var_names to use for computing the hierarchical clustering.
        If `var_names` is given, then `use_rep` and `n_pcs` is ignored.
    use_raw
        Only when `var_names` is not None.
        Use `raw` attribute of `adata` if present.
    cor_method
        correlation method to use.
        Options are 'pearson', 'kendall', and 'spearman'
    linkage_method
        linkage method to use. See :func:`scipy.cluster.hierarchy.linkage`
        for more information.
    optimal_ordering
        Same as the optimal_ordering argument of :func:`scipy.cluster.hierarchy.linkage`
        which reorders the linkage matrix so that the distance between successive
        leaves is minimal.
    key_added
        By default, the dendrogram information is added to
        `.uns[f'dendrogram_{{groupby}}']`.
        Notice that the `groupby` information is added to the dendrogram.
    inplace
        If `True`, adds dendrogram information to `adata.uns[key_added]`,
        else this function returns the information.

    Returns
    -------
    Returns `None` if `inplace=True`, else returns a `dict` with dendrogram information. Sets the following field if `inplace=True`:

    `adata.uns[f'dendrogram_{{group_by}}' | key_added]` : :class:`dict`
        Dendrogram information.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, groupby='bulk_labels')
    >>> sc.pl.dendrogram(adata, groupby='bulk_labels')
    <Axes: >
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    """
    axis, axis_name = _resolve_axis(axis)

    if groupby is not None:
        if isinstance(groupby, str):
            # if not a list, turn into a list
            groupby = [groupby]
        rep_df, categories = _dendrogram_grouped(
            adata,
            groupby,
            n_pcs=n_pcs,
            use_rep=use_rep,
            var_names=var_names,
            use_raw=use_raw,
        )
    else:
        rep_df = pd.DataFrame(
            _choose_representation(
                adata if var_names is None else adata[:, var_names],
                use_rep=use_rep,
                n_pcs=n_pcs,
            )
        )
        categories = rep_df.axes[axis]

    import scipy.cluster.hierarchy as sch
    from scipy.spatial import distance

    corr_matrix = (rep_df.T if axis == 0 else rep_df).corr(method=cor_method)
    corr_condensed = distance.squareform(1 - corr_matrix)
    z_var = sch.linkage(
        corr_condensed, method=linkage_method, optimal_ordering=optimal_ordering
    )
    dendro_info = sch.dendrogram(z_var, labels=list(categories), no_plot=True)

    dat = dict(
        linkage=z_var,
        groupby=groupby,
        use_rep=use_rep,
        cor_method=cor_method,
        linkage_method=linkage_method,
        categories_ordered=dendro_info["ivl"],
        categories_idx_ordered=dendro_info["leaves"],
        dendrogram_info=dendro_info,
        correlation_matrix=corr_matrix.values,
    )

    if not inplace:
        return dat
    key_added = _get_dendrogram_key(key_added, groupby, axis_name=axis_name)
    logg.info(f"Storing dendrogram info using `.uns[{key_added!r}]`")
    adata.uns[key_added] = dat
    return None


def _dendrogram_grouped(
    adata: AnnData,
    groupby: Sequence[str],
    *,
    n_pcs: int | None,
    use_rep: str | None,
    var_names: Sequence[str] | None,
    use_raw: bool | None,
) -> tuple[pd.DataFrame, pd.Index[str]]:
    for group in groupby:
        if group not in adata.obs_keys():
            raise ValueError(
                "groupby has to be a valid observation. "
                f"Given value: {group}, valid observations: {adata.obs_keys()}"
            )
        if not isinstance(adata.obs[group].dtype, CategoricalDtype):
            raise ValueError(
                "groupby has to be a categorical observation. "
                f"Given value: {group}, Column type: {adata.obs[group].dtype}"
            )

    if var_names is None:
        rep_df = pd.DataFrame(
            _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)
        )
        if len(groupby) == 1:
            categorical = adata.obs[groupby[0]]
        else:
            categorical = adata.obs[groupby].apply("_".join, axis=1).astype("category")
        categorical.name = "_".join(groupby)

        rep_df.set_index(categorical, inplace=True)
        categories = rep_df.index.categories
    else:
        from scanpy.plotting._anndata import _prepare_dataframe

        categories, rep_df = _prepare_dataframe(
            adata, var_names, groupby, use_raw=use_raw
        )

        # aggregate values within categories using 'mean'
    rep_df = rep_df.groupby(level=0, observed=True).mean()
    return rep_df, categories


def _get_dendrogram_key(
    dendrogram_key: bool | str | None,
    groupby: str | Sequence[str] | None,
    *,
    axis_name: Literal["obs", "var"],
    adata: AnnData | None = None,
) -> str:
    # the `dendrogram_key` can be a bool an NoneType or the name of the
    # dendrogram key. By default the name of the dendrogram key is 'dendrogram'
    if not isinstance(dendrogram_key, str):
        if groupby is None:
            dendrogram_key = f"dendrogram_{axis_name}"
        elif isinstance(groupby, str):
            dendrogram_key = f"dendrogram_{groupby}"
        elif isinstance(groupby, Sequence):
            dendrogram_key = f'dendrogram_{"_".join(groupby)}'
        else:
            msg = "Either `groupby` or `dendrogram_key` must be specified."
            raise TypeError(msg)

    if adata is None:
        return dendrogram_key

    if dendrogram_key not in adata.uns:
        logg.warning(
            f"dendrogram data not found (using key={dendrogram_key}). "
            "Running `sc.tl.dendrogram` with default parameters. For fine "
            "tuning it is recommended to run `sc.tl.dendrogram` independently."
        )
        dendrogram(adata, groupby, key_added=dendrogram_key)

    if "dendrogram_info" not in adata.uns[dendrogram_key]:
        raise ValueError(
            f"The given dendrogram key ({dendrogram_key!r}) does not contain "
            "valid dendrogram information."
        )

    return dendrogram_key
