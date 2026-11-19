"""Benchmark the aggregation-based plots in Scanpy.

API documentation: <https://scanpy.scverse.org/page/api/plotting.html>.
"""

from __future__ import annotations

import warnings

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy import sparse

import scanpy as sc

mpl.use("Agg")

N_OBS = 100_000
N_VARS = 2_000
N_GROUPS = 15


def _sparse_adata() -> ad.AnnData:
    rng = np.random.default_rng(0)
    x = sparse.random(
        N_OBS, N_VARS, density=0.03, format="csr", dtype=np.float32, random_state=rng
    )
    obs = pd.DataFrame(
        {"group": pd.Categorical(rng.integers(0, N_GROUPS, N_OBS).astype(str))},
        index=[str(i) for i in range(N_OBS)],
    )
    return ad.AnnData(
        x, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(N_VARS)])
    )


class PlottingAggregationSuite:
    """`dotplot` and `matrixplot` should not need a dense cell-level frame."""

    params: tuple[list[int]] = ([20, 200],)
    param_names = ("n_genes",)

    def setup(self, n_genes: int) -> None:
        self.adata = _sparse_adata()
        self.var_names = list(self.adata.var_names[:n_genes])
        # Compile scanpy's aggregation kernels outside of the measured region
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for plot in (sc.pl.dotplot, sc.pl.matrixplot):
                plot(
                    self.adata, self.var_names[:2], "group", return_fig=True, show=False
                )

    def time_dotplot(self, *_) -> None:
        sc.pl.dotplot(self.adata, self.var_names, "group", return_fig=True, show=False)

    def peakmem_dotplot(self, *_) -> None:
        sc.pl.dotplot(self.adata, self.var_names, "group", return_fig=True, show=False)

    def time_matrixplot(self, *_) -> None:
        sc.pl.matrixplot(
            self.adata, self.var_names, "group", return_fig=True, show=False
        )

    def peakmem_matrixplot(self, *_) -> None:
        sc.pl.matrixplot(
            self.adata, self.var_names, "group", return_fig=True, show=False
        )
