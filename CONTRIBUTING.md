# Contributing

Contributions to Scanpy are highly welcome!

## Before filing an issue

* Search the repository (also google) to see if someone has already reported the same issue. This allows contributors to spend less time responding to issues, and more time on adding features!
* Please provide a minimal complete verifiable example for any bug. If you're not sure what this means, check out [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports) by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
* Let us know a bit about your environment. This can be as easy as pasting the results of `sc.logging.print_versions()`.

## Contributing code

## Coding style

We stick to [PEP 8](https://www.python.org/dev/peps/pep-0008) and to this [editorconfig](https://github.com/theislab/scanpy/blob/master/.editorconfig) and *try* to stick to 80-character lines. In some cases, wider lines might improve readability, in most cases, not. Docstrings should always be 80 characters.

### Docs

We use the numpydoc style for writing docstrings. Either take a look at any Scanpy or Numpy function or [here](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

The `Params` abbreviation is a legit replacement for `Parameters`.

The `Returns` section deserves special attention: you can either use the standard numpydoc way of populating it, e.g. as in [`pp.calculate_qc_metrics`](https://scanpy.readthedocs.io/en/latest/api/scanpy.pp.calculate_qc_metrics.html)
```
Returns
-------
type of return value 1
    description of return value 1
type of return value 2
    description of return value 2
```
or you can write plain RST-formatted prose, for instance from [`pp.normalize_total`](https://scanpy.readthedocs.io/en/latest/api/scanpy.pp.normalize_total.html)
```
Returns
-------
Returns dictionary with normalized copies of `adata.X` and `adata.layers`
or updates `adata` with normalized version of the original
`adata.X` and `adata.layers`, depending on `inplace`.
```
or from [`tl.leiden`](https://scanpy.readthedocs.io/en/latest/api/scanpy.tl.leiden.html)
```
Returns
-------
* `adata.obs[key_added]`: Array of dim (number of samples) that stores the subgroup id (`'0'`, `'1'`, ...) for each cell.
* `adata.uns['leiden']['params']`: A dict with the values for the parameters `resolution`, `random_state`, and `n_iterations`.
```
or from [`tl.dpt`](https://scanpy.readthedocs.io/en/latest/api/scanpy.tl.dpt.html)
```
Returns
-------
Depending on `copy`, returns or updates `adata` with the following fields.

If `n_branchings==0`, no field `dpt_groups` will be written.

dpt_pseudotime : `pd.Series` (`adata.obs`, dtype `float`)
    Array of dim (number of samples) that stores the pseudotime of each
    cell, that is, the DPT distance with respect to the root cell.
dpt_groups : `pd.Series` (`adata.obs`, dtype `category`)
    Array of dim (number of samples) that stores the subgroup id ('0',
    '1', ...) for each cell. The groups  typically correspond to
    'progenitor cells', 'undecided cells' or 'branches' of a process.
```

Whether a section is interpreted as `numpydoc` or prose depends on whether the second line of the section is indented or not.


### Tests

Write tests for your functions! See [here](https://github.com/theislab/scanpy/tree/master/scanpy/tests) for examples.
