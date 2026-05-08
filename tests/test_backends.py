from __future__ import annotations

from copy import copy
from inspect import Parameter, signature

import numpy as np
import pytest
from anndata import AnnData

import scanpy as sc
from scanpy._backends import dispatcher
from scanpy._backends import settings as backend_settings
from scanpy.testing import validate_backend


class FakeRapidsBackend:
    name = "rapids_singlecell"
    aliases = ("cuda", "rapids", "rapids-singlecell")

    def normalize_total(
        self,
        adata: AnnData,
        *,
        target_sum: float | None = None,
        fake_param: str | None = None,
    ) -> None:
        """Fake normalize_total backend implementation.

        Parameters
        ----------
        adata
            Annotated data matrix.
        target_sum
            Target total counts.
        fake_param
            Backend-only parameter used by tests.
        """
        sc.pp.normalize_total(adata, target_sum=target_sum, backend="cpu")
        adata.uns["fake_backend_called"] = {
            "function": "normalize_total",
            "fake_param": fake_param,
            "target_sum": target_sum,
        }


DISPATCHED_FUNCTIONS = [
    sc.pp.calculate_qc_metrics,
    sc.pp.filter_cells,
    sc.pp.filter_genes,
    sc.pp.log1p,
    sc.pp.highly_variable_genes,
    sc.pp.normalize_total,
    sc.pp.regress_out,
    sc.pp.scale,
    sc.pp.pca,
    sc.pp.harmony_integrate,
    sc.pp.scrublet,
    sc.pp.scrublet_simulate_doublets,
    sc.pp.neighbors,
    sc.tl.umap,
    sc.tl.tsne,
    sc.tl.diffmap,
    sc.tl.draw_graph,
    sc.tl.embedding_density,
    sc.tl.louvain,
    sc.tl.leiden,
    sc.tl.score_genes,
    sc.tl.score_genes_cell_cycle,
    sc.tl.rank_genes_groups,
    sc.get.aggregate,
]


@pytest.fixture
def fake_rapids_backend():
    registry = dispatcher._registry
    dispatch_impl = dispatcher._dispatch_impl
    old_backend = backend_settings.backend
    old_state = {
        "_backends": copy(registry._backends),
        "_alias_map": copy(registry._alias_map),
        "_load_errors": copy(registry._load_errors),
        "_registration_errors": copy(registry._registration_errors),
        "_warned_untrusted": copy(registry._warned_untrusted),
        "_discovered": registry._discovered,
        "_sig_cache": copy(dispatch_impl._sig_cache),
    }

    backend_settings._backend_var.set("cpu")
    registry._backends.clear()
    registry._alias_map.clear()
    registry._load_errors.clear()
    registry._registration_errors.clear()
    registry._warned_untrusted.clear()
    registry._discovered = True
    registry._register_backend(
        FakeRapidsBackend(),
        entrypoint_name="rapids_singlecell",
        distribution_name="rapids-singlecell",
        object_ref="rapids_singlecell.backends.scanpy:ScanpyBackend",
    )
    dispatch_impl._sig_cache.clear()
    dispatch_impl._update_signatures()

    yield

    backend_settings._backend_var.set(old_backend)
    registry._backends.clear()
    registry._backends.update(old_state["_backends"])
    registry._alias_map.clear()
    registry._alias_map.update(old_state["_alias_map"])
    registry._load_errors.clear()
    registry._load_errors.update(old_state["_load_errors"])
    registry._registration_errors.clear()
    registry._registration_errors.update(old_state["_registration_errors"])
    registry._warned_untrusted.clear()
    registry._warned_untrusted.update(old_state["_warned_untrusted"])
    registry._discovered = old_state["_discovered"]
    dispatch_impl._sig_cache.clear()
    dispatch_impl._sig_cache.update(old_state["_sig_cache"])
    dispatch_impl._update_signatures()


@pytest.mark.parametrize("func", DISPATCHED_FUNCTIONS)
def test_dispatched_functions_have_backend_keyword(func):
    backend = signature(func).parameters["backend"]

    assert backend.kind is Parameter.KEYWORD_ONLY
    assert backend.default is None


def test_settings_resolve_rapids_alias(fake_rapids_backend):
    sc.settings.backend = "cuda"

    assert sc.settings.backend == "rapids_singlecell"
    assert sc.settings.get_backend("rapids") is not None
    assert sc.settings.available_backends() == ["rapids_singlecell"]


def test_use_backend_dispatches_and_restores(fake_rapids_backend):
    adata = AnnData(np.ones((2, 2), dtype=np.float32))

    with sc.settings.use_backend("rapids"):
        assert sc.settings.backend == "rapids_singlecell"
        sc.pp.normalize_total(adata, target_sum=4, fake_param="from-backend")

    assert sc.settings.backend == "cpu"
    assert adata.uns["fake_backend_called"] == {
        "function": "normalize_total",
        "fake_param": "from-backend",
        "target_sum": 4,
    }


def test_call_backend_overrides_settings(fake_rapids_backend):
    adata = AnnData(np.array([[1, 1], [2, 2]], dtype=np.float32))

    with sc.settings.use_backend("cuda"):
        sc.pp.normalize_total(adata, target_sum=1, backend="cpu")

    assert "fake_backend_called" not in adata.uns
    np.testing.assert_allclose(adata.X.sum(axis=1), [1, 1])


def test_backend_only_parameters_are_injected(fake_rapids_backend):
    fake_param = signature(sc.pp.normalize_total).parameters["fake_param"]

    assert fake_param.kind is Parameter.KEYWORD_ONLY
    assert fake_param.default is None


def test_backend_conformance_harness(fake_rapids_backend):
    assert validate_backend("cuda", functions=["normalize_total"]) == {
        "normalize_total": "PASSED",
    }


def test_reserved_gpu_backend_name():
    with pytest.raises(ValueError, match="reserved by scanpy"):
        sc.settings.backend = "gpu"
