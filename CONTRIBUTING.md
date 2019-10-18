Contributing
============

Contributions to Scanpy are highly welcome!

Before filing an issue
----------------------
* Search the repository (also google) to see if someone has already reported the same issue.
  This allows contributors to spend less time responding to issues, and more time adding new features!
* Please provide a minimal complete verifiable example for any bug.
  If you're not sure what this means, check out
  [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports)
  by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
* Let us know about your environment. Environment information is available via: `sc.logging.print_versions()`.

Contributing code
-----------------

We love code contributions! We have a couple guidelines we'd like you to follow though:

### Tests

Please write tests! You can refer to the [existing test suite](https://github.com/theislab/scanpy/tree/master/scanpy/tests) for examples. If you haven't written tests before, Software Carpentry has an [in-depth guide](http://katyhuff.github.io/python-testing/) on the topic.

Test are run by issuing the command `pytest` from the root of the repository. `pytest` as well as a few other testing dependencies can be installed by running `pip install ".[test]"` from the repository root, or `pip install scanpy[test]`.

### Coding style
New code should follow [Black][] and Scanpy’s [EditorConfig][],
so using an editor/IDE with support for both is helpful.

[Black]: https://black.readthedocs.io/en/stable/the_black_code_style.html
[EditorConfig]: https://github.com/theislab/scanpy/blob/master/.editorconfig

### Docs and type annotations
We use the numpydoc style for writing docstrings.
Look at [`sc.tl.louvain`][] as an example for everything mentioned here:

The `Params` abbreviation is a legit replacement for `Parameters`.

To document parameter types use type annotations on function parameters.
Use the [`typing`][] module for containers, e.g. `Sequence`s (like `list`),
`Iterable`s (like `set`), and `Mapping`s (like `dict`). Always specify
what these contain, e.g. `{'a': (1, 2)}` → `Mapping[str, Tuple[int, int]]`.
If you can’t use one of those, use a concrete class like `AnnData`.
If your parameter only accepts an enumeration of strings, specify them like so:
`Literal['elem-1', 'elem-2']`.

The `Returns` section deserves special attention:
There are three types of return sections – prose, tuple, and a mix of both.

1. Prose is for simple cases.
2. Tuple return sections are formatted like parameters.
   Other than in numpydoc, each tuple is first characterized by the identifier
   and *not* by its type. Provide type annotation in the function header.
3. Mix of prose and tuple is relevant in complicated cases,
   e.g. when you want to describe that you
   *added something as annotation to an `AnnData` object*.

[`sc.tl.louvain`]: https://github.com/theislab/scanpy/blob/a811fee0ef44fcaecbde0cad6336336bce649484/scanpy/tools/_louvain.py#L22-L90
[`typing`]: https://docs.python.org/3/library/typing.html

#### Examples
For simple cases, use prose as in [`pp.normalize_total`][].

```rst
Returns
-------
Returns dictionary with normalized copies of ``adata.X`` and ``adata.layers``
or updates ``adata`` with normalized versions of the original
``adata.X`` and ``adata.layers``, depending on ``inplace``.
```

You can use the standard numpydoc way of populating it,
e.g. as in [`pp.calculate_qc_metrics`][].
If you use a plain type name here, a link will be created.

```rst
Returns
-------
one_identifier : some_module.some_type
    Description.
second_identifier : another.module.and_type
    Description 2.
```

Many functions also just modify parts of the passed AnnData object,
like e.g. [`tl.dpt`][].
You can then combine prose and lists to best describe what happens.

```rst
Returns
-------
Depending on `copy`, returns or updates `adata` with the following fields.

If `n_branchings==0`, no field `dpt_groups` will be written.

dpt_pseudotime : :class:`~pandas.Series` (``adata.obs``, dtype ``float``)
    Array of dim (number of samples) that stores the pseudotime of each
    cell, that is, the DPT distance with respect to the root cell.
dpt_groups : :class:`pandas.Series` (``adata.obs``, dtype ``category``)
    Array of dim (number of samples) that stores the subgroup id ('0',
    '1', ...) for each cell. The groups  typically correspond to
    'progenitor cells', 'undecided cells' or 'branches' of a process.
```

[`pp.normalize_total`]: https://scanpy.readthedocs.io/en/latest/api/scanpy.pp.normalize_total.html
[`pp.calculate_qc_metrics`]: https://scanpy.readthedocs.io/en/latest/api/scanpy.pp.calculate_qc_metrics.html
[`tl.dpt`]: https://scanpy.readthedocs.io/en/latest/api/scanpy.tl.dpt.html

### Performance

We defer loading a few modules until they’re first needed.
If you want realistic performance measures,
be sure to import them before running scanpy functions:

- Check the list in `test_deferred_imports()` from [`scanpy.tests.test_performance`][]
- Everything in [`scanpy.external`][] wraps a 3rd party import.

[`scanpy.tests.test_performance`]: https://github.com/theislab/scanpy/blob/master/scanpy/tests/test_performance.py
[`scanpy.external`]: https://scanpy.readthedocs.io/en/stable/external/
