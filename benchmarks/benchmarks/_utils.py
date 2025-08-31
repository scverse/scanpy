from __future__ import annotations

import itertools
import warnings
from functools import cache
from typing import TYPE_CHECKING

import numpy as np
import pooch
from anndata import concat
from asv_runner.benchmarks.mark import skip_for_params

import scanpy as sc
from scanpy._compat import CSRBase

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from collections.abc import Set as AbstractSet
    from typing import Literal, Protocol, TypeVar

    from anndata import AnnData

    from scanpy._compat import CSCBase

    C = TypeVar("C", bound=Callable)

    class ParamSkipper(Protocol):
        def __call__(self, **skipped: AbstractSet) -> Callable[[C], C]: ...

    Dataset = Literal["pbmc68k_reduced", "pbmc3k", "bmmc", "lung93k"]
    KeyX = Literal[None, "off-axis"]
    KeyCount = Literal["counts", "counts-off-axis"]


@cache
def _pbmc68k_reduced() -> AnnData:
    """A small datasets with a dense `.X`."""  # noqa: D401
    adata = sc.datasets.pbmc68k_reduced()
    assert isinstance(adata.X, np.ndarray)
    assert not np.isfortran(adata.X)

    # raw has the same number of genes, so we can use it for counts
    # it doesn’t actually contain counts for some reason, but close enough
    assert isinstance(adata.raw.X, CSRBase)
    adata.layers["counts"] = adata.raw.X.toarray(order="C")
    mapper = dict(
        percent_mito="pct_counts_mt",
        n_counts="total_counts",
    )
    adata.obs.rename(columns=mapper, inplace=True)
    return adata


def pbmc68k_reduced() -> AnnData:
    return _pbmc68k_reduced().copy()


@cache
def _pbmc3k() -> AnnData:
    adata = sc.datasets.pbmc3k()
    assert isinstance(adata.X, CSRBase)
    adata.layers["counts"] = adata.X.astype(np.int32, copy=True)
    sc.pp.log1p(adata)
    return adata


def pbmc3k() -> AnnData:
    return _pbmc3k().copy()


@cache
def _bmmc(n_obs: int = 4000) -> AnnData:
    registry = pooch.create(
        path=pooch.os_cache("pooch"),
        base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    registry.load_registry_from_doi()
    samples = {smp: f"{smp}_filtered_feature_bc_matrix.h5" for smp in ("s1d1", "s1d3")}
    adatas = {}

    for sample_id, filename in samples.items():
        path = registry.fetch(filename)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"Variable names are not unique")
            sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        sc.pp.subsample(sample_adata, n_obs=n_obs // len(samples))
        adatas[sample_id] = sample_adata

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Observation names are not unique")
        adata = concat(adatas, label="sample")
    adata.obs_names_make_unique()

    assert isinstance(adata.X, CSRBase)
    adata.layers["counts"] = adata.X.astype(np.int32, copy=True)
    sc.pp.log1p(adata)
    adata.obs["n_counts"] = adata.layers["counts"].sum(axis=1).A1
    return adata


def bmmc(n_obs: int = 400) -> AnnData:
    return _bmmc(n_obs).copy()


@cache
def _lung93k() -> AnnData:
    path = pooch.retrieve(
        url="https://figshare.com/ndownloader/files/45788454",
        known_hash="md5:4f28af5ff226052443e7e0b39f3f9212",
    )
    adata = sc.read_h5ad(path)
    assert isinstance(adata.X, CSRBase)
    adata.layers["counts"] = adata.X.astype(np.int32, copy=True)
    sc.pp.log1p(adata)
    return adata


def lung93k() -> AnnData:
    return _lung93k().copy()


def to_off_axis(x: np.ndarray | CSRBase) -> np.ndarray | CSCBase:
    if isinstance(x, CSRBase):
        return x.tocsc()
    if isinstance(x, np.ndarray):
        assert not np.isfortran(x)
        return x.copy(order="F")
    msg = f"Unexpected type {type(x)}"
    raise TypeError(msg)


def _get_dataset_raw(dataset: Dataset) -> tuple[AnnData, str | None]:
    match dataset:
        case "pbmc68k_reduced":
            adata, batch_key = pbmc68k_reduced(), None
        case "pbmc3k":
            adata, batch_key = pbmc3k(), None  # can’t use this with batches
        case "bmmc":
            # TODO: allow specifying bigger variant
            adata, batch_key = bmmc(400), "sample"
        case "lung93k":
            adata, batch_key = lung93k(), "PatientNumber"
        case _:
            msg = f"Unknown dataset {dataset}"
            raise AssertionError(msg)

    # add off-axis layers
    adata.layers["off-axis"] = to_off_axis(adata.X)
    adata.layers["counts-off-axis"] = to_off_axis(adata.layers["counts"])

    # add mitochondrial gene and pre-compute qc metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    assert adata.var["mt"].sum() > 0, "no MT genes in dataset"
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    return adata, batch_key


def get_dataset(dataset: Dataset, *, layer: KeyX = None) -> tuple[AnnData, str | None]:
    adata, batch_key = _get_dataset_raw(dataset)
    if layer is not None:
        adata.X = adata.layers.pop(layer)
    return adata, batch_key


def get_count_dataset(
    dataset: Dataset, *, layer: KeyCount = "counts"
) -> tuple[AnnData, str | None]:
    adata, batch_key = _get_dataset_raw(dataset)

    adata.X = adata.layers.pop(layer)
    # remove indicators that X was transformed
    adata.uns.pop("log1p", None)

    return adata, batch_key


def param_skipper(
    param_names: Sequence[str], params: tuple[Sequence[object], ...]
) -> ParamSkipper:
    """Create a decorator that will skip all combinations that contain any of the given parameters.

    Examples
    --------
    >>> param_names = ["letters", "numbers"]
    >>> params = [["a", "b"], [3, 4, 5]]
    >>> skip_when = param_skipper(param_names, params)

    >>> @skip_when(letters={"a"}, numbers={3})
    ... def func(a, b):
    ...     print(a, b)
    >>> run_as_asv_benchmark(func)
    b 4
    b 5

    """

    def skip(**skipped: AbstractSet) -> Callable[[C], C]:
        skipped_combs = [
            tuple(record.values())
            for record in (
                dict(zip(param_names, vals, strict=True))
                for vals in itertools.product(*params)
            )
            if any(v in skipped.get(n, set()) for n, v in record.items())
        ]
        # print(skipped_combs, file=sys.stderr)
        return skip_for_params(skipped_combs)

    return skip
