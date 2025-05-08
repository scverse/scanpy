"""Computes a dendrogram based on a given categorical observation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from pandas.api.types import CategoricalDtype

from .. import logging as logg
from .._compat import old_positionals
from .._utils import _doc_params, raise_not_implemented_error_if_backed_type
from ..neighbors._doc import doc_n_pcs, doc_use_rep
from ._utils import _choose_representation

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

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
def dendrogram(  # noqa: PLR0913
    adata: AnnData,
    groupby: str | Sequence[str],
    *,
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
    """Compute a hierarchical clustering for the given `groupby` categories.

    By default, the PCA representation is used unless `.X`
    has less than 50 variables.

    Alternatively, a list of `var_names` (e.g. genes) can be given.

    Average values of either `var_names` or components are used
    to compute a correlation matrix.

    The hierarchical clustering can be visualized using
    :func:`scanpy.pl.dendrogram` or multiple other visualizations
    that can include a dendrogram: :func:`~scanpy.pl.matrixplot`,
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
    {n_pcs}
    {use_rep}
    var_names
        List of var_names to use for computing the hierarchical clustering.
        If `var_names` is given, then `use_rep` and `n_pcs` are ignored.
    use_raw
        Only when `var_names` is not None.
        Use `raw` attribute of `adata` if present.
    cor_method
        Correlation method to use.
        Options are 'pearson', 'kendall', and 'spearman'
    linkage_method
        Linkage method to use. See :func:`scipy.cluster.hierarchy.linkage`
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
    >>> sc.tl.dendrogram(adata, groupby="bulk_labels")
    >>> sc.pl.dendrogram(adata, groupby="bulk_labels")  # doctest: +SKIP
    <Axes: >
    >>> markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
    >>> sc.pl.dotplot(adata, markers, groupby="bulk_labels", dendrogram=True)

    """
    raise_not_implemented_error_if_backed_type(adata.X, "dendrogram")
    if isinstance(groupby, str):
        # if not a list, turn into a list
        groupby = [groupby]
    for group in groupby:
        if group not in adata.obs_keys():
            msg = (
                "groupby has to be a valid observation. "
                f"Given value: {group}, valid observations: {adata.obs_keys()}"
            )
            raise ValueError(msg)
        if not isinstance(adata.obs[group].dtype, CategoricalDtype):
            msg = (
                "groupby has to be a categorical observation. "
                f"Given value: {group}, Column type: {adata.obs[group].dtype}"
            )
            raise ValueError(msg)

    if var_names is None:
        rep_df = pd.DataFrame(
            _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)
        )
        categorical = adata.obs[groupby[0]]
        if len(groupby) > 1:
            for group in groupby[1:]:
                # create new category by merging the given groupby categories
                categorical = (
                    categorical.astype(str) + "_" + adata.obs[group].astype(str)
                ).astype("category")
        categorical.name = "_".join(groupby)

        rep_df.set_index(categorical, inplace=True)
        categories: pd.Index = rep_df.index.categories
    else:
        gene_names = adata.raw.var_names if use_raw else adata.var_names
        from ..plotting._anndata import _prepare_dataframe

        categories, rep_df = _prepare_dataframe(
            adata, gene_names, groupby, use_raw=use_raw
        )

    # aggregate values within categories using 'mean'
    mean_df = (
        rep_df.groupby(level=0, observed=True)
        .mean()
        .loc[categories]  # Fixed ordering for pandas < 2
    )

    import scipy.cluster.hierarchy as sch
    from scipy.spatial import distance

    corr_matrix = mean_df.T.corr(method=cor_method).clip(-1, 1)
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

    if inplace:
        if key_added is None:
            key_added = f"dendrogram_{'_'.join(groupby)}"
        logg.info(f"Storing dendrogram info using `.uns[{key_added!r}]`")
        adata.uns[key_added] = dat
    else:
        return dat
