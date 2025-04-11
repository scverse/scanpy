from __future__ import annotations

import importlib
import os
from collections import defaultdict
from inspect import Parameter, signature
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

import pytest
from anndata import AnnData

# CLI is locally not imported by default but on travis it is?
import scanpy.cli
from scanpy._utils import _import_name, descend_classes_and_funcs

if TYPE_CHECKING:
    from types import FunctionType
    from typing import Any

mod_dir = Path(scanpy.__file__).parent
proj_dir = mod_dir.parent


api_module_names = [
    "sc",
    "sc.pp",
    "sc.tl",
    "sc.pl",
    "sc.experimental.pp",
    "sc.external.pp",
    "sc.external.tl",
    "sc.external.pl",
    "sc.external.exporting",
    "sc.get",
    "sc.logging",
    # "sc.neighbors",  # Not documented
    "sc.datasets",
    "sc.queries",
    "sc.metrics",
]
api_modules = {
    mod_name: _import_name(f"scanpy{mod_name.removeprefix('sc')}")
    for mod_name in api_module_names
}


# get all exported functions that aren’t re-exports from anndata
api_functions = [
    pytest.param(func, f"{mod_name}.{name}", id=f"{mod_name}.{name}")
    for mod_name, mod in api_modules.items()
    for name in sorted(mod.__all__)
    if callable(func := getattr(mod, name)) and func.__module__.startswith("scanpy.")
]


@pytest.fixture
def in_project_dir():
    wd_orig = Path.cwd()
    os.chdir(proj_dir)
    try:
        yield proj_dir
    finally:
        os.chdir(wd_orig)


@pytest.mark.xfail(reason="TODO: unclear if we want this to totally match, let’s see")
def test_descend_classes_and_funcs():
    funcs = set(descend_classes_and_funcs(scanpy, "scanpy"))
    assert {p.values[0] for p in api_functions} == funcs


@pytest.mark.filterwarnings("error::FutureWarning:.*Import anndata.*")
def test_import_future_anndata_import_warning():
    import scanpy

    importlib.reload(scanpy)


def param_is_pos(p: Parameter) -> bool:
    return p.kind in {
        Parameter.POSITIONAL_ONLY,
        Parameter.POSITIONAL_OR_KEYWORD,
    }


def is_deprecated(f: FunctionType) -> bool:
    # TODO: use deprecated decorator instead
    # https://github.com/scverse/scanpy/issues/2505
    return f.__name__ in {
        "normalize_per_cell",
        "filter_genes_dispersion",
    }


class ExpectedSig(TypedDict):
    first_name: str
    copy_default: Any
    return_ann: str | None


copy_sigs: defaultdict[str, ExpectedSig | None] = defaultdict(
    lambda: ExpectedSig(first_name="adata", copy_default=False, return_ann=None)
)
# full exceptions
copy_sigs["sc.external.tl.phenograph"] = None  # external
copy_sigs["sc.pp.filter_genes_dispersion"] = None  # deprecated
copy_sigs["sc.pp.filter_cells"] = None  # unclear `inplace` situation
copy_sigs["sc.pp.filter_genes"] = None  # unclear `inplace` situation
copy_sigs["sc.pp.subsample"] = None  # returns indices along matrix
copy_sigs["sc.pp.sample"] = None  # returns indices along matrix
# partial exceptions: “data” instead of “adata”
copy_sigs["sc.pp.log1p"]["first_name"] = "data"
copy_sigs["sc.pp.normalize_per_cell"]["first_name"] = "data"
copy_sigs["sc.pp.pca"]["first_name"] = "data"
copy_sigs["sc.pp.scale"]["first_name"] = "data"
copy_sigs["sc.pp.sqrt"]["first_name"] = "data"
# other partial exceptions
copy_sigs["sc.pp.normalize_total"]["return_ann"] = copy_sigs[
    "sc.experimental.pp.normalize_pearson_residuals"
]["return_ann"] = "AnnData | dict[str, np.ndarray] | None"
copy_sigs["sc.external.pp.magic"]["copy_default"] = None


@pytest.mark.parametrize(("f", "qualname"), api_functions)
def test_sig_conventions(f, qualname):
    sig = signature(f)

    # TODO: replace the following check with lint rule for all funtions eventually
    if not is_deprecated(f):
        n_pos = sum(1 for p in sig.parameters.values() if param_is_pos(p))
        assert n_pos <= 3, "Public functions should have <= 3 positional parameters"

    first_param = next(iter(sig.parameters.values()), None)
    if first_param is None:
        return

    if first_param.name == "adata":
        assert first_param.annotation in {"AnnData", AnnData}
    elif first_param.name == "data":
        assert first_param.annotation.startswith("AnnData |")
    elif first_param.name in {"filename", "path"}:
        assert first_param.annotation == "Path | str"

    # Test if functions with `copy` follow conventions
    if (copy_param := sig.parameters.get("copy")) is not None and (
        expected_sig := copy_sigs[qualname]
    ) is not None:
        s = ExpectedSig(
            first_name=first_param.name,
            copy_default=copy_param.default,
            return_ann=sig.return_annotation,
        )
        expected_sig = expected_sig.copy()
        if expected_sig["return_ann"] is None:
            expected_sig["return_ann"] = f"{first_param.annotation} | None"
        assert s == expected_sig
        if not is_deprecated(f):
            assert not param_is_pos(copy_param)


def getsourcefile(obj):
    """inspect.getsourcefile, but supports singledispatch."""
    from inspect import getsourcefile

    if wrapped := getattr(obj, "__wrapped__", None):
        return getsourcefile(wrapped)

    return getsourcefile(obj)


def getsourcelines(obj):
    """inspect.getsourcelines, but supports singledispatch."""
    from inspect import getsourcelines

    if wrapped := getattr(obj, "__wrapped__", None):
        return getsourcelines(wrapped)

    return getsourcelines(obj)
