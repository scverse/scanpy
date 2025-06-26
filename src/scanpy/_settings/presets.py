from __future__ import annotations

import inspect
import re
from contextlib import contextmanager
from enum import StrEnum, auto
from functools import cached_property, partial, wraps
from typing import TYPE_CHECKING, Literal, NamedTuple, TypeVar

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Mapping


NT = TypeVar("NT", bound=NamedTuple)

__all__ = [
    "DETest",
    "HVGFlavor",
    "HVGPreset",
    "LeidenFlavor",
    "LeidenPreset",
    "PcaPreset",
    "Preset",
    "RankGenesGroupsPreset",
]


DETest = Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]
HVGFlavor = Literal["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"]
LeidenFlavor = Literal["leidenalg", "igraph"]


class HVGPreset(NamedTuple):
    flavor: HVGFlavor
    return_df: bool


class PcaPreset(NamedTuple):
    key_added: str | None


class RankGenesGroupsPreset(NamedTuple):
    method: DETest
    mask_var: str | None


class LeidenPreset(NamedTuple):
    flavor: LeidenFlavor


preset_postprocessors: list[Callable[[], None]] = []


def named_tuple_non_defaults(
    nt: NamedTuple,
) -> Generator[tuple[str, object], None, None]:
    cls = type(nt)
    for param in cls._fields:
        value = getattr(nt, param)
        if param not in cls._field_defaults or value != cls._field_defaults[param]:
            yield param, value


def postprocess_preset_prop(
    prop: cached_property[NT], get_map: Callable[[], Mapping[Preset, NT]]
) -> None:
    map = get_map()

    map_type = inspect.signature(get_map).return_annotation
    m = re.fullmatch(r"Mapping\[Preset, (.*)\]", map_type)
    assert m is not None
    value_type = m[1]

    added_doc = "\n".join(
        ":attr:`{name}`\n    Defaults: {defaults}".format(
            name=k.name,
            defaults=", ".join(
                f"`{param}={default!r}`"
                for param, default in named_tuple_non_defaults(params)
            )
            or "none",
        )
        for k, params in map.items()
    )

    prop.__doc__ = f"{prop.__doc__}\n\n{added_doc}"
    prop.func.__annotations__["return"] = value_type


def preset_property(get_map: Callable[[], Mapping[Preset, NT]]) -> cached_property[NT]:
    @wraps(get_map)
    def get(self: Preset) -> NT:
        return get_map()[self]

    prop = cached_property(get)
    preset_postprocessors.append(partial(postprocess_preset_prop, prop, get_map))
    return prop


class Preset(StrEnum):
    """Presets for :attr:`scanpy.settings.preset`.

    See properties below for details.
    """

    ScanpyV1 = auto()
    """Scanpy 1.*’s default settings."""

    ScanpyV2Preview = auto()
    """Scanpy 2.*’s feature default settings. (Preview: subject to change!)"""

    SeuratV5 = auto()
    """Try to match Seurat 5.* as closely as possible."""

    @preset_property
    def highly_variable_genes() -> Mapping[Preset, HVGPreset]:
        """Flavor for :func:`~scanpy.pp.highly_variable_genes`."""
        return {
            Preset.ScanpyV1: HVGPreset(flavor="seurat", return_df=False),
            Preset.ScanpyV2Preview: HVGPreset(flavor="seurat_v3_paper", return_df=True),
            Preset.SeuratV5: HVGPreset(flavor="seurat_v3_paper", return_df=True),
        }

    @preset_property
    def pca() -> Mapping[Preset, PcaPreset]:
        """Settings for :func:`~scanpy.pp.pca`."""  # noqa: D401
        return {
            Preset.ScanpyV1: PcaPreset(key_added=None),
            Preset.ScanpyV2Preview: PcaPreset(key_added="pca"),
            Preset.SeuratV5: PcaPreset(key_added="pca"),
        }

    @preset_property
    def rank_genes_groups() -> Mapping[Preset, RankGenesGroupsPreset]:
        """Correlation method for :func:`~scanpy.tl.rank_genes_groups`."""
        return {
            Preset.ScanpyV1: RankGenesGroupsPreset(method="t-test", mask_var=None),
            Preset.ScanpyV2Preview: RankGenesGroupsPreset(
                method="wilcoxon", mask_var=None
            ),
            Preset.SeuratV5: RankGenesGroupsPreset(
                method="wilcoxon", mask_var="highly_variable"
            ),
        }

    @preset_property
    def leiden() -> Mapping[Preset, LeidenPreset]:
        """Flavor for :func:`~scanpy.tl.leiden`."""
        return {
            Preset.ScanpyV1: LeidenPreset(flavor="leidenalg"),
            Preset.ScanpyV2Preview: LeidenPreset(flavor="igraph"),
            Preset.SeuratV5: LeidenPreset(flavor="leidenalg"),
        }

    @contextmanager
    def override(self, preset: Preset) -> Generator[Preset, None, None]:
        """Temporarily override :attr:`scanpy.settings.preset`.

        >>> import scanpy as sc
        >>> sc.settings.preset = sc.Preset.ScanpyV1
        >>> with sc.settings.preset.override(sc.Preset.SeuratV5):
        ...     sc.settings.preset
        <Preset.SeuratV5: 'seuratv5'>
        >>> sc.settings.preset
        <Preset.ScanpyV1: 'scanpyv1'>
        """
        from scanpy import settings

        settings.preset = preset
        try:
            yield self
        finally:
            settings.preset = self


for postprocess in preset_postprocessors:
    postprocess()
