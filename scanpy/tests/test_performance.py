def test_deferred_imports(imported_modules):
    slow_to_import = {
        'umap',  # neighbors, tl.umap
        'seaborn',  # plotting
        'sklearn.metrics',  # neighbors
        'scipy.stats',  # tools._embedding_density
        'networkx',  # diffmap, paga, plotting._utils
        # TODO: 'matplotlib.pyplot',
        # TODO (maybe): 'numba',
    }
    falsely_imported = slow_to_import & imported_modules
    assert not falsely_imported
