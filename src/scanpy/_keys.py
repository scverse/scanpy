from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, assert_never, cast, overload

from ._compat import warn
from ._settings import Default, Preset, settings

if TYPE_CHECKING:
    from anndata import AnnData

    from ._settings.presets import BasicEmbeddingPreset


__all__ = [
    "_Embedding",
    "_EmbeddingKeys",
    "_PcaKeys",
    "_embedding_keys",
    "_existing_preset_keys",
]


type _Embedding = Literal["pca", "tsne", "umap", "draw_graph", "diffmap"]


@dataclass(frozen=True)
class _EmbeddingKeys:
    uns: str
    obsm: str


@dataclass(frozen=True)
class _PcaKeys(_EmbeddingKeys):
    varm: str


@overload
def _embedding_keys(
    embedding: Literal["pca"],
    key_added: str | None | Default | Preset = ...,
) -> _PcaKeys: ...
@overload
def _embedding_keys(
    embedding: Literal["tsne", "umap", "diffmap"],
    key_added: str | None | Default | Preset = ...,
) -> _EmbeddingKeys: ...
@overload
def _embedding_keys(
    embedding: Literal["draw_graph"],
    key_added: str | None | Default | Preset = ...,
    *,
    layout: str,
    key_added_ext: str | None = ...,
) -> _EmbeddingKeys: ...
def _embedding_keys(  # noqa: PLR0911
    embedding: _Embedding,
    key_added: str | None | Default | Preset = Default(),
    *,
    layout: str = "",
    key_added_ext: str | None = None,
) -> _EmbeddingKeys:
    if isinstance(key_added, Default):
        key_added = settings.preset
    if isinstance(key_added, Preset):
        key_added = cast(
            "BasicEmbeddingPreset", getattr(key_added, embedding)
        ).key_added
    match embedding, key_added:
        case "draw_graph", _:
            return _draw_graph_keys(
                key_added, layout=layout, key_added_ext=key_added_ext
            )
        case "pca", str():
            return _PcaKeys(key_added, key_added, key_added)
        case _, str():
            return _EmbeddingKeys(key_added, key_added)
        case "pca", None:
            return _PcaKeys("pca", "X_pca", "PCs")
        case "tsne", None:
            return _EmbeddingKeys("tsne", "X_tsne")
        case "umap", None:
            return _EmbeddingKeys("umap", "X_umap")
        case "diffmap", None:
            return _EmbeddingKeys("diffmap_evals", "X_diffmap")
        case _:
            assert_never(embedding)


def _draw_graph_keys(
    key_added: str | None, *, layout: str, key_added_ext: str | None
) -> _EmbeddingKeys:
    if key_added_ext is not None:
        msg = "Passing `key_added_ext` is deprecated, use `key_added`'s template functionality instead."
        warn(msg, category=FutureWarning)
        suffix = key_added_ext
    else:
        suffix = layout
    if key_added is None:
        return _EmbeddingKeys("draw_graph", f"X_draw_graph_{suffix}")
    formatted = key_added.format(layout=suffix)
    return _EmbeddingKeys(formatted, formatted)


@overload
def _existing_preset_keys(
    adata: AnnData, embedding: Literal["pca"]
) -> _PcaKeys | None: ...
@overload
def _existing_preset_keys(
    adata: AnnData, embedding: Literal["tsne", "umap", "diffmap"]
) -> _EmbeddingKeys | None: ...
@overload
def _existing_preset_keys(
    adata: AnnData, embedding: Literal["draw_graph"], *, layout: str
) -> _EmbeddingKeys | None: ...
def _existing_preset_keys(
    adata: AnnData,
    embedding: _Embedding,
    *,
    layout: str = "",
) -> _EmbeddingKeys | None:
    for preset in (Preset.ScanpyV1, Preset.ScanpyV2Preview):
        if embedding == "draw_graph":
            result = _embedding_keys("draw_graph", preset, layout=layout)
        else:
            result = _embedding_keys(embedding, preset)  # type: ignore[call-overload]
        if result.obsm in adata.obsm:
            return result
    return None
