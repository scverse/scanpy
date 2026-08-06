"""Backend dispatch integration for optional Scanpy accelerators."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scverse_backends import BackendDispatcher

if TYPE_CHECKING:
    from collections.abc import Callable

dispatcher = BackendDispatcher(
    entrypoint_group="scanpy.backends",
    host_name="scanpy",
    trusted_backends={
        "rapids-singlecell": {
            "aliases": ["cuda", "rapids", "rapids_singlecell"],
            "distributions": [
                "rapids-singlecell",
                "rapids-singlecell-cu12",
                "rapids-singlecell-cu13",
            ],
            "package": "rapids-singlecell",
        },
    },
    reserved_backends={
        "gpu": (
            "Use 'cuda' for the RAPIDS backend. The generic 'gpu' selector is "
            "reserved so future GPU backends can coexist without ambiguity."
        ),
    },
)

dispatched_functions: list[Callable[..., object]] = []


def backend_dispatch[**P, R](func: Callable[P, R]) -> Callable[P, R]:
    """Decorate a public Scanpy function for backend dispatch."""
    dispatched = dispatcher.backend_dispatch(func)
    dispatched_functions.append(dispatched)
    return dispatched


settings = dispatcher.settings
get_backend = dispatcher.get_backend
available_backend_names = dispatcher.available_backend_names
discover = dispatcher.discover

__all__ = [
    "available_backend_names",
    "backend_dispatch",
    "discover",
    "dispatched_functions",
    "dispatcher",
    "get_backend",
    "settings",
]
