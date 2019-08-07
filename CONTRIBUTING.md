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
We stick to [PEP 8](https://www.python.org/dev/peps/pep-0008) and this
[editorconfig](https://github.com/theislab/scanpy/blob/master/.editorconfig)
and *try* to stick to 80-character lines.
In some cases, wider lines might improve readability, in most cases, not.
Docstrings should always be 80 characters.

### Docs
We use the numpydoc style for writing docstrings.
Either take a look at any Scanpy or Numpy function or
[here](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

The `Params` abbreviation is a legit replacement for `Parameters`.

The `Returns` section deserves special attention:
There are three types of return sections – prose, tuple, and a mix of both.

1. Prose is for simple cases.
2. Tuple return sections are formatted like parameters.
   Other than in numpydoc, each tuple is first characterized by the identifier name
   and *not* by its type. You can provide type annotation in the function header
   or by separation with a colon, as in parameters.
3. Mix of prose and tuple is relevant in complicated cases,
   e.g. when you want to describe that you *added something as annotation to an `AnnData` object*.

#### Examples
For simple cases, use prose as in [`pp.normalize_total`](https://scanpy.readthedocs.io/en/latest/api/scanpy.pp.normalize_total.html)

```rst
Returns
-------
Returns dictionary with normalized copies of ``adata.X`` and ``adata.layers``
or updates ``adata`` with normalized versions of the original
``adata.X`` and ``adata.layers``, depending on ``inplace``.
```

You can use the standard numpydoc way of populating it, e.g. as in
[`pp.calculate_qc_metrics`](https://scanpy.readthedocs.io/en/latest/api/scanpy.pp.calculate_qc_metrics.html).
If you just use a plain type name here, there will be an automatically created link.

```rst
Returns
-------
one_identifier : some_module.some_type
    Description.
second_identifier : another.module.and_type
    Description 2.
```

Many functions also just modify parts of the passed AnnData object,
like e.g. [`tl.dpt`](https://scanpy.readthedocs.io/en/latest/api/scanpy.tl.dpt.html).
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

### Performance

We defer loading a few modules until they’re first needed.
If you want realistic performance measures, be sure to import them before running scanpy functions:

- Check the list in `test_deferred_imports()` from [`scanpy/tests/test_performance.py`](https://github.com/theislab/scanpy/blob/master/scanpy/tests/test_performance.py)
- Everything in [`scanpy.external`](https://scanpy.readthedocs.io/en/stable/external/) wraps a 3rd party import.
