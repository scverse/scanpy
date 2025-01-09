# Documentation

(building-the-docs)=

## Building the docs

To build the docs, run `hatch run docs:build`.
Afterwards, you can run `hatch run docs:open` to open {file}`docs/_build/html/index.html`.

Your browser and Sphinx cache docs which have been built previously.
Sometimes these caches are not invalidated when you've updated the docs.
If docs are not updating the way you expect, first try "force reloading" your browser page – e.g. reload the page without using the cache.
Next, if problems persist, clear the sphinx cache (`hatch run docs:clean`) and try building them again.

(adding-to-the-docs)=

## Adding to the docs

For any user-visible changes, please make sure a note has been added to the release notes using [`hatch run towncrier:create`][towncrier create].
When asked for “Issue number (`+` if none)”, enter the *PR number* instead.

Once you've added a new function to the documentation, you'll need to make sure there is a link somewhere in the documentation site pointing to it.
This should be added to `docs/api.md` under a relevant heading.

For tutorials and more in depth examples, consider adding a notebook to the [scanpy-tutorials][] repository.

The tutorials are tied to this repository via a submodule.
To update the submodule, run `git submodule update --remote` from the root of the repository.
Subsequently, commit and push the changes in a PR.
This should be done before each release to ensure the tutorials are up to date.

[towncrier create]: https://towncrier.readthedocs.io/en/stable/tutorial.html#creating-news-fragments
[scanpy-tutorials]: https://github.com/scverse/scanpy-tutorials/

## docstrings format

We use the numpydoc style for writing docstrings.
We'd primarily suggest looking at existing docstrings for examples, but the [napolean guide to numpy style docstrings][] is also a great source.
If you're unfamiliar with the reStructuredText (rST) markup format, check out the [Sphinx rST primer][].

Some key points:

- We have some custom sphinx extensions activated. When in doubt, try to copy the style of existing docstrings.
- We autopopulate type information in docstrings when possible, so just add the type information to signatures.
- When docs exist in the same file as code, line length restrictions still apply. In files which are just docs, go with a sentence per line (for easier `git diff`s).
- Check that the docs look like what you expect them too! It's easy to forget to add a reference to function, be sure it got added and looks right.

Look at [sc.tl.louvain](https://github.com/scverse/scanpy/blob/a811fee0ef44fcaecbde0cad6336336bce649484/scanpy/tools/_louvain.py#L22-L90) as an example for everything mentioned here.

[napolean guide to numpy style docstrings]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy
[sphinx rst primer]: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

### Plots in docstrings

One of the most useful things you can include in a docstring is examples of how the function should be used.
These are a great way to demonstrate intended usage and give users a template they can copy and modify.
We're able to include the plots produced by these snippets in the rendered docs using [matplotlib's plot directive][].
For examples of this, see the `Examples` sections of {func}`~scanpy.pl.dotplot` or {func}`~scanpy.pp.calculate_qc_metrics`.

Note that anything in these sections will need to be run when the docs are built, so please keep them computationally light.

- If you need computed features (e.g. an embedding, differential expression results) load data that has this precomputed.
- Try to re-use datasets, this reduces the amount of data that needs to be downloaded to the CI server.

[matplotlib's plot directive]: https://matplotlib.org/devel/plot_directive.html

### `Params` section

The `Params` abbreviation is a legit replacement for `Parameters`.

To document parameter types use type annotations on function parameters.
These will automatically populate the docstrings on import, and when the documentation is built.

Use the python standard library types (defined in {mod}`collections.abc` and {mod}`typing` modules) for containers, e.g.
{class}`~collections.abc.Sequence`s (like `list`),
{class}`~collections.abc.Iterable`s (like `set`), and
{class}`~collections.abc.Mapping`s (like `dict`).
Always specify what these contain, e.g. `{'a': (1, 2)}` → `Mapping[str, Tuple[int, int]]`.
If you can’t use one of those, use a concrete class like `AnnData`.
If your parameter only accepts an enumeration of strings, specify them like so: `Literal['elem-1', 'elem-2']`.

### `Returns` section

There are three types of return sections – prose, tuple, and a mix of both.

1. Prose is for simple cases.
2. Tuple return sections are formatted like parameters. Other than in numpydoc, each tuple is first characterized by the identifier and *not* by its type. Provide type annotation in the function header.
3. Mix of prose and tuple is relevant in complicated cases, e.g. when you want to describe that you *added something as annotation to an \`AnnData\` object*.

#### Examples

For simple cases, use prose as in {func}`~scanpy.pp.normalize_total`:

```rst
Returns
-------
Returns dictionary with normalized copies of `adata.X` and `adata.layers`
or updates `adata` with normalized versions of the original
`adata.X` and `adata.layers`, depending on `inplace`.
```

For tuple return values, you can use the standard numpydoc way of populating it,
e.g. as in {func}`~scanpy.pp.calculate_qc_metrics`.
Do not add types in the docstring, but specify them in the function signature:

```python
def myfunc(...) -> tuple[int, str]:
    """
    ...
    Returns
    -------
    one_identifier
        Description.
    second_identifier
        Description 2.
    """
    ...
```

Many functions also just modify parts of the passed AnnData object, like e.g. {func}`~scanpy.tl.dpt`.
You can then combine prose and lists to best describe what happens:

```rst
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
```
