from __future__ import annotations

import enum
import inspect
import re
from contextlib import contextmanager
from dataclasses import dataclass
from functools import cached_property, partial, wraps
from importlib.metadata import distributions, requires
from typing import TYPE_CHECKING, Literal, NamedTuple

from packaging.requirements import Requirement

from .._utils._doctests import doctest_needs

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Mapping


__all__ = [
    "DETest",
    "Default",
    "HVGFlavor",
    "HVGPreset",
    "LeidenFlavor",
    "LeidenPreset",
    "PcaPreset",
    "Preset",
    "RankGenesGroupsPreset",
]


type DETest = Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]
type HVGFlavor = Literal["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"]
type LeidenFlavor = Literal["leidenalg", "igraph"]


@dataclass
class Default:
    repr: str | None = None
    preset: tuple[str, str] | None = None

    def __post_init__(self):
        if self.preset and self.repr:
            msg = "Cannot provide both preset and repr."
            raise TypeError(msg)

    def __repr__(self) -> str:
        if not self.preset:
            return self.repr or "default"
        import scanpy as sc

        value = self._get_value(sc.settings.preset)
        suffix = (
            " – changes in 2.0"
            if sc.settings.preset is Preset.ScanpyV1
            and value != self._get_value(Preset.ScanpyV2Preview)
            else ""
        )
        return f"{value!r} (sc.settings.preset={str(sc.settings.preset)!r}{suffix})"

    def _get_value(self, preset: Preset) -> object:
        if not self.preset:
            msg = "preset is not set"
            raise AssertionError(msg)
        func, param = self.preset
        params = getattr(preset, func)
        return getattr(params, param)


class HVGPreset(NamedTuple):
    flavor: HVGFlavor
    return_df: bool


class PcaPreset(NamedTuple):
    key_added: str | None


class RankGenesGroupsPreset(NamedTuple):
    method: DETest
    mask_var: str | None


class ScalePreset(NamedTuple):
    zero_center: bool | None


class ScoreGenesPreset(NamedTuple):
    ctrl_as_ref: bool


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


def postprocess_preset_prop[NT: NamedTuple](
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


def preset_property[NT: NamedTuple](
    get_map: Callable[[], Mapping[Preset, NT]],
) -> cached_property[NT]:
    @wraps(get_map)
    def get(self: Preset) -> NT:
        return get_map()[self]

    prop = cached_property(get)
    preset_postprocessors.append(partial(postprocess_preset_prop, prop, get_map))  # noqa: F821
    return prop


class Preset(enum.StrEnum):
    """Presets for :attr:`scanpy.settings.preset`.

    See properties below for details.
    """

    @staticmethod
    def _generate_next_value_(
        name: str, start: int, count: int, last_values: list[str]
    ):
        # lower-kebap-case
        return "-".join(part.lower() for part in re.split(r"(?=[A-Z])", name) if part)

    ScanpyV1 = enum.auto()
    """: Scanpy 1.*’s default settings."""

    ScanpyV2Preview = enum.auto()
    """: Scanpy 2.*’s feature default settings. (Preview: subject to change!)"""

    @preset_property
    def highly_variable_genes() -> Mapping[Preset, HVGPreset]:
        """Flavor for :func:`~scanpy.pp.highly_variable_genes`."""
        return {
            Preset.ScanpyV1: HVGPreset(flavor="seurat", return_df=False),
            Preset.ScanpyV2Preview: HVGPreset(flavor="seurat_v3_paper", return_df=True),
        }

    @preset_property
    def pca() -> Mapping[Preset, PcaPreset]:
        """Settings for :func:`~scanpy.pp.pca`."""  # noqa: D401
        return {
            Preset.ScanpyV1: PcaPreset(key_added=None),
            Preset.ScanpyV2Preview: PcaPreset(key_added="pca"),
        }

    @preset_property
    def rank_genes_groups() -> Mapping[Preset, RankGenesGroupsPreset]:
        """Correlation method for :func:`~scanpy.tl.rank_genes_groups`."""
        return {
            Preset.ScanpyV1: RankGenesGroupsPreset(method="t-test", mask_var=None),
            Preset.ScanpyV2Preview: RankGenesGroupsPreset(
                method="wilcoxon", mask_var=None
            ),
        }

    @preset_property
    def scale() -> Mapping[Preset, ScalePreset]:
        """Flavor for :func:`~scanpy.pp.scale`."""
        return {
            Preset.ScanpyV1: ScalePreset(zero_center=True),
            Preset.ScanpyV2Preview: ScalePreset(zero_center=None),
        }

    @preset_property
    def score_genes() -> Mapping[Preset, ScoreGenesPreset]:
        """Settings for :func:`~scanpy.tl.score_genes`."""  # noqa: D401
        return {
            Preset.ScanpyV1: ScoreGenesPreset(ctrl_as_ref=True),
            Preset.ScanpyV2Preview: ScoreGenesPreset(ctrl_as_ref=False),
        }

    @preset_property
    def leiden() -> Mapping[Preset, LeidenPreset]:
        """Flavor for :func:`~scanpy.tl.leiden`."""
        return {
            Preset.ScanpyV1: LeidenPreset(flavor="leidenalg"),
            Preset.ScanpyV2Preview: LeidenPreset(flavor="igraph"),
        }

    @contextmanager
    @doctest_needs("igraph")
    @doctest_needs("scikit-misc")
    def override(self, preset: Preset) -> Generator[Preset, None, None]:
        """Temporarily override :attr:`scanpy.settings.preset`.

        >>> import scanpy as sc
        >>> sc.settings.preset = sc.Preset.ScanpyV1
        >>> with sc.settings.preset.override(sc.Preset.ScanpyV2Preview):
        ...     sc.settings.preset
        <Preset.ScanpyV2Preview: 'scanpy-v2-preview'>
        >>> sc.settings.preset
        <Preset.ScanpyV1: 'scanpy-v1'>
        """
        from scanpy import settings

        settings.preset = preset
        try:
            yield self
        finally:
            settings.preset = self

    def check(self) -> None:
        """Check if requirements for preset are met."""
        match self:
            case self.ScanpyV1:
                return
            case self.ScanpyV2Preview:
                if not (missing := _missing_scanpy2_deps()):
                    return
                missing_str = ", ".join(f"‘{m.name}’" for m in missing)
                msg = (
                    f"Setting preset to {Preset.ScanpyV2Preview!r} requires optional "
                    f"dependencies that are not installed: {missing_str}. "
                    "Install them with: pip install `scanpy[scanpy2]`"
                )
                raise ImportError(msg)


def _missing_scanpy2_deps() -> list[Requirement]:
    dist_names = {d.name for d in distributions()}
    return [
        r
        for r in map(Requirement, requires("scanpy") or ())
        if r.marker
        and r.marker.evaluate({"extra": "scanpy2"})
        and r.name not in dist_names
    ]


for postprocess in preset_postprocessors:
    postprocess()
del postprocess, preset_postprocessors
