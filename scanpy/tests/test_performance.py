import json
import sys
from subprocess import run


def descend(profimp_data, modules, path):
    module = profimp_data["module"]
    path = [*path, module]
    if module in modules:
        yield " → ".join(e for e in path if e is not None)
        modules.remove(module)
    for child in profimp_data["children"]:
        yield from descend(child, modules, path)


def get_import_paths(modules):
    proc = run(
        [sys.executable, "-m", "profimp.main", "import scanpy"],
        capture_output=True,
        check=True,
    )
    data = json.loads(proc.stdout)
    return descend(data, set(modules), [])


def test_deferred_imports(imported_modules):
    slow_to_import = {
        'umap',  # neighbors, tl.umap
        'seaborn',  # plotting
        'sklearn.metrics',  # neighbors
        'networkx',  # diffmap, paga, plotting._utils
        # TODO: 'matplotlib.pyplot',
        # TODO (maybe): 'numba',
    }
    falsely_imported = slow_to_import & imported_modules

    assert not falsely_imported, "\n".join(get_import_paths(falsely_imported))
