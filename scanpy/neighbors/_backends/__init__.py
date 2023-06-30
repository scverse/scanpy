from __future__ import annotations

from collections import ChainMap
from typing import Union

from . import umap, rapids, third_party
from ._common import mappings, select_backend, Transformer

# TODO gauss
_Backend = Union[umap._Backend, rapids._Backend, third_party._Backend]
_Algorithm = Union[umap._Algorithm, rapids._Algorithm, third_party._Algorithm]

_backends: ChainMap[_Backend, frozenset[_Algorithm]] = ChainMap(
    umap.BACKENDS, rapids.BACKENDS, third_party.BACKENDS
)
BACKENDS, ALGORITHMS = mappings(_backends)


def get_transformer(
    algorithm: _Algorithm, backend: _Backend | None = None
) -> type[Transformer]:
    """
    Get the transformer for a given algorithm and backend.

    Infers the backend when that is unambiguously possible.
    """
    backend = select_backend(ALGORITHMS, algorithm, backend)
    if backend in umap.BACKENDS:
        return umap.UMAPKNNTransformer
    if backend in rapids.BACKENDS:
        return rapids.RapidsKNNTransformer
    return third_party.get_transformer(algorithm, backend)
