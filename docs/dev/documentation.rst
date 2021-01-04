Documentation
=============

.. _building-the-docs:

Building the docs
-----------------

Dependencies for building the documentation for scanpy can be installed with `pip install -e "scanpy[doc]"`

To build the docs, enter the `docs` directory and run `make html`. After this process completes you can take a look at the docs by opening `scanpy/docs/_build/html/index.html`.

Your browser and Sphinx cache docs which have been built previously.
Sometimes these caches are not invalidated when you've updated the docs.
If docs are not updating the way you expect, first try "force reloading" your browser page – e.g. reload the page without using the cache.
Next, if problems persist, clear the sphinx cache and try building them again (`make clean` from `docs` directory).


Adding to the docs
------------------

For any user-visible changes, please make sure a note has been added to `docs/release-latest.rst` so we can credit you!
We recommend waiting on this until your PR is close to done since this section often causes merge conflicts.

Once you've added a new function to the documentation, you'll need to make sure there is a link somewhere in the documentation site pointing to it.
For computational methods, this should be added to `docs/api/index.rst` under a relevant heading.
For plotting functions, add these to the module docstring of the plotting module at `scanpy/pl/__init__.py`.

For tutorials and more in depth examples, consider adding a notebook to `scanpy-tutorials <https://github.com/theislab/scanpy-tutorials/>`__.

docstrings format
-----------------

We use the numpydoc style for writing docstrings.
We'd primarily suggest looking at existing docstrings for examples, but the `napolean guide to numpy style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy>`__ is also a great source.
If you're unfamiliar with the reStructuredText (`rst`) markup format, `Sphinx has a useful primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__.

Some key points:

* We have some custom sphinx extensions activated. When in doubt, try to copy the style of existing docstrings.
* We autopopulate type information in docstrings when possible, so just add the type information to signatures.
* When docs exist in the same file as code, line length restrictions still apply. In files which are just docs, go with a sentence per line (for easier `git diff`\ s).
* Check that the docs look like what you expect them too! It's easy to forget to add a reference to function, be sure it got added and looks right.

Look at `sc.tl.louvain <https://github.com/theislab/scanpy/blob/a811fee0ef44fcaecbde0cad6336336bce649484/scanpy/tools/_louvain.py#L22-L90>`__ as an example for everything mentioned here:

`Params` section
~~~~~~~~~~~~~~~~

The `Params` abbreviation is a legit replacement for `Parameters`.

To document parameter types use type annotations on function parameters.
These will automatically populate the docstrings on import, and when the documentation is built.

Use the python standard library types (defined in `collections.abc <https://docs.python.org/3/library/collections.abc.html>`__ and `typing <https://docs.python.org/3/library/typing.html>`__ modules) for containers, e.g. `Sequence`\ s (like `list`), `Iterable`\ s (like `set`), and `Mapping`\ s (like `dict`).
Always specify what these contain, e.g. `{'a': (1, 2)}` → `Mapping[str, Tuple[int, int]]`.
If you can’t use one of those, use a concrete class like `AnnData`.
If your parameter only accepts an enumeration of strings, specify them like so: `Literal['elem-1', 'elem-2']`.

`Returns` section
~~~~~~~~~~~~~~~~~

There are three types of return sections – prose, tuple, and a mix of both.

1. Prose is for simple cases.
2. Tuple return sections are formatted like parameters. Other than in numpydoc, each tuple is first characterized by the identifier and *not* by its type. Provide type annotation in the function header.
3. Mix of prose and tuple is relevant in complicated cases, e.g. when you want to describe that you *added something as annotation to an `AnnData` object*.

Examples
^^^^^^^^

For simple cases, use prose as in
:func:`~scanpy.pp.normalize_total`

.. code:: rst

   Returns
   -------
   Returns dictionary with normalized copies of `adata.X` and `adata.layers`
   or updates `adata` with normalized versions of the original
   `adata.X` and `adata.layers`, depending on `inplace`.

You can use the standard numpydoc way of populating it, e.g. as in
:func:`~scanpy.pp.calculate_qc_metrics`.
If you use a plain type name here, a link will be created.

.. code:: rst

   Returns
   -------
   one_identifier : some_module.some_type
       Description.
   second_identifier : another.module.and_type
       Description 2.

Many functions also just modify parts of the passed AnnData object, like e.g. :func:`~scanpy.tl.dpt`.
You can then combine prose and lists to best describe what happens.

.. code:: rst

   Returns
   -------
   Depending on `copy`, returns or updates `adata` with the following fields.

   If `n_branchings==0`, no field `dpt_groups` will be written.

   dpt_pseudotime : :class:`~pandas.Series` (`adata.obs`, dtype `float`)
       Array of dim (number of samples) that stores the pseudotime of each
       cell, that is, the DPT distance with respect to the root cell.
   dpt_groups : :class:`pandas.Series` (`adata.obs`, dtype `category`)
       Array of dim (number of samples) that stores the subgroup id ('0',
       '1', ...) for each cell. The groups  typically correspond to
       'progenitor cells', 'undecided cells' or 'branches' of a process.
