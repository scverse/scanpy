from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from numpy.typing import NDArray

    from .._compat import CSBase


def rename_groups(
    adata: AnnData,
    restrict_key: str,
    *,
    key_added: str | None,
    restrict_categories: Iterable[str],
    restrict_indices: NDArray[np.bool_],
    groups: NDArray,
) -> pd.Series[str]:
    key_added = f"{restrict_key}_R" if key_added is None else key_added
    all_groups = adata.obs[restrict_key].astype("U")
    prefix = "-".join(restrict_categories) + ","
    new_groups = [prefix + g for g in groups.astype("U")]
    all_groups.iloc[restrict_indices] = new_groups
    return all_groups


def restrict_adjacency(
    adata: AnnData,
    restrict_key: str,
    *,
    restrict_categories: Sequence[str],
    adjacency: CSBase,
) -> tuple[CSBase, NDArray[np.bool_]]:
    if not isinstance(restrict_categories[0], str):
        msg = "You need to use strings to label categories, e.g. '1' instead of 1."
        raise ValueError(msg)
    for c in restrict_categories:
        if c not in adata.obs[restrict_key].cat.categories:
            msg = f"{c!r} is not a valid category for {restrict_key!r}"
            raise ValueError(msg)
    restrict_indices = adata.obs[restrict_key].isin(restrict_categories).values
    adjacency = adjacency[restrict_indices, :]
    adjacency = adjacency[:, restrict_indices]
    return adjacency, restrict_indices
