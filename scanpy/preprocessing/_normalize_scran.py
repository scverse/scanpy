from typing import Collection, Optional, Union

import numpy as np
from anndata import AnnData
from anndata.utils import Index

from .. import logging as logg


def normalize_scran(
    adata: AnnData,
    *,
    sizes: Collection[int] = np.arange(21, 102, 5),
    clusters: Union[str, Collection[Union[str, int]], None] = None,
    ref_clust: Optional[Union[str, int]] = None,
    max_cluster_size: int = 3000,
    # positive: bool = True,
    scaling: Union[str, Collection[float], None] = None,
    min_mean: Union[int, float] = 1,
    subset_gene: Union[str, Index, None] = None,
    inplace: bool = True,
):
    """\
    One sentence to describe it.

    More details, ...

    Parameters
    ----------
    adata
        ...
    sizes
        ...
    ...
    positive
        Should linear inverse models be used to enforce positive estimates?
    ...

    Returns
    -------
    """
    if isinstance(clusters, str):
        clusters = adata.obs[clusters]
    if isinstance(scaling, str):
        scaling = adata.obs[scaling]
    if isinstance(subset_gene, str):
        subset_gene = adata.var[subset_gene]

    ...

    return adata if not inplace else None
