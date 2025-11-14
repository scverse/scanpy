"""Add `array-support` directive."""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from itertools import product
from typing import TYPE_CHECKING, overload

if TYPE_CHECKING:
    from collections.abc import Collection, Generator
    from typing import Literal


__all__ = ["ArrayType", "DaskArray", "Numpy", "ScipySparse", "parse"]


class ArrayType(ABC):
    def rst(self) -> str:
        return f":class:`{self}`"

    @abstractmethod
    def __hash__(self) -> int: ...


@dataclass(unsafe_hash=True, frozen=True)
class Numpy(ArrayType):
    def __str__(self) -> str:
        return "numpy.ndarray"

    def rst(self) -> str:
        return f":class:`{self}`"


@dataclass(unsafe_hash=True, frozen=True)
class ScipySparse(ArrayType):
    format: Literal["csr", "csc"]
    container: Literal["array", "matrix"]

    def __str__(self) -> str:
        return f"scipy.sparse.{self.format}_{self.container}"


type Inner = Numpy | ScipySparse


@dataclass(unsafe_hash=True, frozen=True)
class DaskArray(ArrayType):
    chunk: Inner

    def __str__(self) -> str:
        return f"dask.array.Array[{self.chunk}]"

    def rst(self) -> str:
        return rf":class:`dask.array.Array`\ \[{self.chunk.rst()}\]"


@overload
def parse(
    include: Collection[str],
    exclude: Collection[str] = (),
    *,
    inner: Literal[False] = False,
) -> Generator[ArrayType]: ...
@overload
def parse(
    include: Collection[str], exclude: Collection[str] = (), *, inner: Literal[True]
) -> Generator[Inner]: ...
def parse(
    include: Collection[str], exclude: Collection[str] = (), *, inner: bool = False
) -> Generator[ArrayType]:
    if exclude:
        excluded = dict.fromkeys(parse(exclude)).keys()
        yield from (t for t in parse(include) if t not in excluded)
        return

    inner_includes = [i for i in include if not i.startswith("da")]
    for t in include:
        if (match := re.fullmatch(r"([^\[]+)(?:\[(.+)\])?", t)) is None:
            msg = f"invalid {t!r}"
            raise ValueError(msg)
        mod, tags = match.groups("")
        if mod == "da" and inner:
            msg = "Can’t nest dask arrays"
            raise ValueError(msg)
        tags = set(re.split(r",(?![^\[]+\])", tags)) if tags else set()
        yield from _parse_mod(mod, tags, inner_includes=inner_includes)


def _parse_mod(
    mod: str, tags: set[str], *, inner_includes: Collection[str]
) -> Generator[ArrayType]:
    match mod:
        case "np":
            if tags:
                msg = f"`np` takes no tags {tags!r}"
                raise ValueError(msg)
            yield Numpy()
        case "sp":
            if tags - {"csr", "csc", "array", "matrix"}:
                msg = f"invalid tags {tags!r}"
                raise ValueError(msg)
            for format, container in product(("csr", "csc"), ("array", "matrix")):
                if tags & {"csr", "csc"} and format not in tags:
                    continue
                if tags & {"array", "matrix"} and container not in tags:
                    continue
                yield ScipySparse(format=format, container=container)
        case "da":
            for chunk in parse(tags if tags else inner_includes, inner=True):
                yield DaskArray(chunk=chunk)
        case _:
            msg = f"invalid module {mod!r}"
            raise ValueError(msg)
