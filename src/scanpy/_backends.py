"""Backend dispatch integration for optional Scanpy accelerators."""

from __future__ import annotations

from scverse_backends import BackendDispatcher

dispatcher = BackendDispatcher(
    entrypoint_group="scanpy.backends",
    host_name="scanpy",
    trusted_backends={
        "rapids_singlecell": {
            "aliases": ["cuda", "rapids", "rapids-singlecell"],
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

backend_dispatch = dispatcher.backend_dispatch
settings = dispatcher.settings
get_backend = dispatcher.get_backend
available_backend_names = dispatcher.available_backend_names
discover = dispatcher.discover

__all__ = [
    "available_backend_names",
    "backend_dispatch",
    "discover",
    "dispatcher",
    "get_backend",
    "settings",
]
