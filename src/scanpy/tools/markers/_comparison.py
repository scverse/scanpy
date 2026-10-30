from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ... import _utils
from ..._compat import CSBase
from ..._utils import raise_not_implemented_error_if_backed_type
from ...get import _check_mask

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import Literal

    from anndata import AnnData
    from numpy.typing import NDArray

    from ..._compat import DaskArray
    from ...get.get import Mask


def _expm1_func(adata: AnnData) -> Callable[[NDArray], NDArray]:
    if (base := adata.uns.get("log1p", {}).get("base")) is None:
        return np.expm1

    def expm1(x: NDArray) -> NDArray:
        return np.expm1(x * np.log(base))

    return expm1


@dataclass(frozen=True, eq=False)
class _Comparison:
    """Which cells are compared, on which matrix."""

    x: NDArray[np.number] | CSBase | DaskArray
    var_names: pd.Index
    labels: pd.Categorical
    groups_order: NDArray
    groups_masks_obs: NDArray[np.bool]
    ireference: int | None
    expm1_func: Callable[[NDArray], NDArray]

    @property
    def n_groups(self) -> int:
        return self.groups_masks_obs.shape[0]

    @property
    def grouping_mask(self) -> NDArray[np.bool]:
        return np.asarray(pd.Series(self.labels).isin(self.groups_order))

    @property
    def grouping(self) -> pd.Categorical:
        return self.labels[self.grouping_mask]

    @property
    def reference(self) -> str:
        if self.ireference is None:
            return "rest"
        return str(self.groups_order[self.ireference])

    @classmethod
    def legacy(
        cls,
        adata: AnnData,
        groups: Iterable[str] | Literal["all"],
        groupby: str,
        *,
        mask_var: NDArray[np.bool] | None = None,
        reference: str = "rest",
        use_raw: bool = True,
        layer: str | None = None,
    ) -> _Comparison:
        groups_order, groups_masks_obs = _utils.select_groups(adata, groups, groupby)
        _check_group_sizes(adata.obs[groupby], groups_order)

        adata_comp = adata
        if layer is not None:
            if use_raw:
                msg = "Cannot specify `layer` and have `use_raw=True`."
                raise ValueError(msg)
            x = adata_comp.layers[layer]
        else:
            if use_raw and adata.raw is not None:
                adata_comp = adata.raw
            x = adata_comp.X
        raise_not_implemented_error_if_backed_type(x, "rank_genes_groups")

        # for correct getnnz calculation
        if isinstance(x, CSBase):
            x.eliminate_zeros()

        var_names = adata_comp.var_names
        if mask_var is not None:
            x = x[:, mask_var]
            var_names = var_names[mask_var]

        return cls(
            x=x,
            var_names=var_names,
            labels=adata.obs[groupby].array,
            groups_order=groups_order,
            groups_masks_obs=groups_masks_obs,
            ireference=_reference_index(groups_order, reference),
            expm1_func=_expm1_func(adata),
        )

    @classmethod
    def from_adata(
        cls,
        adata: AnnData,
        groupby: str,
        *,
        groups: Iterable[str] | Literal["all"] = "all",
        reference: str = "rest",
        layer: str | None = None,
        mask: Mask | None = None,
        func_name: str,
    ) -> _Comparison:
        labels = pd.Categorical(adata.obs[groupby])
        groups_order = _resolve_groups(labels, groups, reference)
        codes = pd.Index(labels.categories).get_indexer(groups_order)
        groups_masks_obs = labels.codes[None, :] == codes[:, None]
        _check_group_sizes(pd.Series(labels), groups_order)

        x = adata.X if layer is None else adata.layers[layer]
        raise_not_implemented_error_if_backed_type(x, func_name)
        var_names = adata.var_names
        if (mask_var := _check_mask(adata, mask, "var")) is not None:
            x = x[:, mask_var]
            var_names = var_names[mask_var]

        return cls(
            x=x,
            var_names=var_names,
            labels=labels,
            groups_order=groups_order,
            groups_masks_obs=groups_masks_obs,
            ireference=_reference_index(groups_order, reference),
            expm1_func=_expm1_func(adata),
        )


def _resolve_groups(
    labels: pd.Categorical, groups: Iterable[str] | Literal["all"], reference: str
) -> NDArray:
    categories = labels.categories
    if reference != "rest" and reference not in categories:
        msg = f"reference = {reference} needs to be one of groupby = {categories.tolist()}."
        raise ValueError(msg)
    if isinstance(groups, str) and groups == "all":
        return categories.to_numpy()
    if isinstance(groups, str | int):
        msg = "Specify a sequence of groups"
        raise ValueError(msg)
    selected = [str(g) if isinstance(g, int) else g for g in groups]
    if missing := [g for g in selected if g not in categories]:
        msg = (
            f"Groups {missing} are not categories of `groupby`: {categories.tolist()}."
        )
        raise ValueError(msg)
    if reference != "rest" and reference not in selected:
        selected.append(reference)
    if len(set(selected) - {reference}) == 0:
        msg = "`groups` needs to contain at least one group other than `reference`."
        raise ValueError(msg)
    return categories[categories.get_indexer(selected)].to_numpy()


def _check_group_sizes(labels: pd.Series, groups_order: NDArray) -> None:
    # Singlet groups cause division by zero errors
    invalid_groups_selected = set(groups_order) & set(
        labels.value_counts().loc[lambda x: x < 2].index
    )
    if len(invalid_groups_selected) > 0:
        msg = (
            f"Could not calculate statistics for groups {', '.join(map(str, invalid_groups_selected))} "
            "since they contain fewer than two cells."
        )
        raise ValueError(msg)


def _reference_index(groups_order: NDArray, reference: str) -> int | None:
    if reference == "rest":
        return None
    return np.where(groups_order == reference)[0][0]
