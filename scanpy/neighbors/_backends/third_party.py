from __future__ import annotations

from importlib import import_module
from types import MappingProxyType, ModuleType
from typing import Literal, Union
from collections.abc import Mapping, Generator

from sklearn.neighbors import KNeighborsTransformer


_Eponymous = Literal['annoy', 'faiss', 'nmslib', 'pynndescent']
_Backend = Union[_Eponymous, Literal['sklearn']]
_Algorithm = Union[_Eponymous, Literal['ball_tree', 'kd_tree', 'brute']]

BACKENDS: Mapping[_Backend, frozenset[_Algorithm]] = MappingProxyType(
    {
        'sklearn': frozenset({'ball_tree', 'kd_tree', 'brute'}),
        'annoy': frozenset({'annoy'}),
        'faiss': frozenset({'faiss'}),
        'nmslib': frozenset({'nmslib'}),
        'pynndescent': frozenset({'pynndescent'}),
    }
)
# TODO: Are there algos supported by multiple backends?
ALGORITHMS: Mapping[_Algorithm, _Backend] = MappingProxyType(
    {algo: backend for backend, algos in BACKENDS.items() for algo in algos}
)


def import_backend(backend: _Backend) -> ModuleType:
    return import_module(f'sklearn_ann.kneighbors.{backend}')


def get_transformer(algorithm: _Algorithm) -> KNeighborsTransformer:
    if (backend := ALGORITHMS.get(algorithm)) is None:
        raise ValueError(f'Unknown algorithm: {algorithm} is not in {set(ALGORITHMS)}')
    mod = import_backend(backend)
    transformer_cls: type[KNeighborsTransformer] = next(
        cls
        for name, cls in vars(mod).items()
        if name.lower().startswith(algorithm.replace('_', ''))
    )
    return transformer_cls()


def available_backends() -> Generator[_Backend, None, None]:
    for backend in BACKENDS:
        try:
            import_backend(backend)
        except ImportError:
            pass
        else:
            yield backend
