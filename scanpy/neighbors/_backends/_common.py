from __future__ import annotations

import operator
from collections.abc import Mapping, Iterable, Collection
from functools import reduce
from typing import Any, TypeVar, Protocol
from types import MappingProxyType

import numpy as np
from scipy.sparse import csr_matrix


BE = TypeVar("BE")
ALG = TypeVar("ALG")


class Transformer(Protocol):
    def __init__(self, *, n_neighbors: int, **kwargs) -> None:
        ...

    def fit(self, X: np.ndarray | csr_matrix, y: Any = None) -> Transformer:
        ...

    def transform(self, X: np.ndarray | csr_matrix | None = None) -> csr_matrix:
        ...

    def fit_transform(self, X: np.ndarray | csr_matrix, y: Any = None) -> csr_matrix:
        ...


class TransformerChecksMixin:
    def _transform_checks(self, X, *fitted_props, **check_params):
        from sklearn.utils.validation import check_is_fitted

        if X is not None:
            X = self._validate_data(X, reset=False, **check_params)
        check_is_fitted(self, *fitted_props)
        return X


def mappings(
    backend2algo: Mapping[BE, Iterable[ALG]]
) -> tuple[Mapping[BE, frozenset[ALG]], Mapping[ALG, frozenset[BE]]]:
    """Freeze backend-to-algorithms mapping and create reverse mapping."""
    be2alg = MappingProxyType(
        {backend: frozenset(algos) for backend, algos in backend2algo.items()}
    )
    all_algos: frozenset[ALG] = reduce(operator.or_, backend2algo.values(), frozenset())
    alg2be = MappingProxyType(
        {
            algo: frozenset(
                backend for backend, algos in backend2algo.items() if algo in algos
            )
            for algo in all_algos
        }
    )
    return be2alg, alg2be


def select_backend(
    algs: Mapping[ALG, Collection[BE]],
    algorithm: ALG,
    backend: BE | None = None,
) -> BE:
    if backend is not None:
        return backend

    if (backends := algs.get(algorithm)) is None:
        msg = f'Unknown algorithm: {algorithm} is not in {set(algs)}'
        raise ValueError(msg)
    if len(backends) > 1:
        msg = f'Algorithm {algorithm} is supported by multiple backends {backends}, specify one'
        raise ValueError(msg)
    return next(iter(backends))
