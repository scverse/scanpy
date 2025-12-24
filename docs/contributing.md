# Contributing guide

This document aims at summarizing the most important information for getting you started on contributing to this project.
We assume that you are already familiar with git and with making pull requests on GitHub.

For more extensive tutorials, that also cover the absolute basics,
please refer to other resources such as the [pyopensci tutorials][],
the [scientific Python tutorials][], or the [scanpy developer guide][].

[pyopensci tutorials]: https://www.pyopensci.org/learn.html
[scientific Python tutorials]: https://learn.scientific-python.org/development/tutorials/
[scanpy developer guide]: https://scanpy.readthedocs.io/en/latest/dev/index.html

:::{tip} The *hatch* project manager

We highly recommend to familiarize yourself with [`hatch`][hatch].
Hatch is a Python project manager that

- manages virtual environments, separately for development, testing and building the documentation.
  Separating the environments is useful to avoid dependency conflicts.
- allows to run tests locally in different environments (e.g. different python versions)
- allows to run tasks defined in `pyproject.toml`, e.g. to build documentation.

While the project is setup with `hatch` in mind,
it is still possible to use different tools to manage dependencies, such as `uv` or `pip`.

:::

[hatch]: https://hatch.pypa.io/latest/

## Installing dev dependencies

In addition to the packages needed to _use_ this package,
you need additional python packages to [run tests](#writing-tests) and [build the documentation](#docs-building).

:::::{tabs}
::::{group-tab} Hatch

On the command line, you typically interact with hatch through its command line interface (CLI).
Running one of the following commands will automatically resolve the environments for testing and
building the documentation in the background:

```bash
hatch test  # defined in the table [tool.hatch.envs.hatch-test] in pyproject.toml
hatch run docs:build  # defined in the table [tool.hatch.envs.docs]
```

When using an IDE such as VS Code,
you’ll have to point the editor at the paths to the virtual environments manually.
The environment you typically want to use as your main development environment is the `hatch-test`
environment with the latest Python version.

To get a list of all environments for your projects, run

```bash
hatch env show -i
```

This will list “Standalone” environments and a table of “Matrix” environments like the following:

```
+------------+---------+--------------------------+----------+---------------------------------+-------------+
| Name       | Type    | Envs                     | Features | Dependencies                    | Scripts     |
+------------+---------+--------------------------+----------+---------------------------------+-------------+
| hatch-test | virtual | hatch-test.py3.11-stable | dev      | coverage-enable-subprocess==1.0 | cov-combine |
|            |         | hatch-test.py3.14-stable | test     | coverage[toml]~=7.4             | cov-report  |
|            |         | hatch-test.py3.14-pre    |          | pytest-mock~=3.12               | run         |
|            |         |                          |          | pytest-randomly~=3.15           | run-cov     |
|            |         |                          |          | pytest-rerunfailures~=14.0      |             |
|            |         |                          |          | pytest-xdist[psutil]~=3.5       |             |
|            |         |                          |          | pytest~=8.1                     |             |
+------------+---------+--------------------------+----------+---------------------------------+-------------+
```

From the `Envs` column, select the environment name you want to use for development.
In this example, it would be `hatch-test.py3.14-stable`.

Next, create the environment with

```bash
hatch env create hatch-test.py3.14-stable
```

Then, obtain the path to the environment using

```bash
hatch env find hatch-test.py3.14-stable
```

In case you are using VScode, now open the command palette (Ctrl+Shift+P) and search for `Python: Select Interpreter`.
Choose `Enter Interpreter Path` and paste the path to the virtual environment from above.

In this future, this may become easier through a hatch vscode extension.

::::

::::{group-tab} uv

A popular choice for managing virtual environments is [uv][].
The main disadvantage compared to hatch is that it supports only a single environment per project at a time,
which requires you to mix the dependencies for running tests and building docs.
This can have undesired side-effects,
such as requiring to install a lower version of a library your project depends on,
only because an outdated sphinx plugin pins an older version.

To initalize a virtual environment in the `.venv` directory of your project, simply run

```bash
uv sync --all-extras
```

The `.venv` directory is typically automatically discovered by IDEs such as VS Code.

::::

::::{group-tab} Pip

Pip is nowadays mostly superseded by environment manager such as [hatch][].
However, for the sake of completeness, and since it’s ubiquitously available,
we describe how you can manage environments manually using `pip`:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,test,doc]"
```

The `.venv` directory is typically automatically discovered by IDEs such as VS Code.

::::
:::::

[hatch environments]: https://hatch.pypa.io/latest/tutorials/environment/basic-usage/
[uv]: https://docs.astral.sh/uv/

## Code-style

This package uses [pre-commit][] to enforce consistent code-styles.
On every commit, pre-commit checks will either automatically fix issues with the code, or raise an error message.

To enable pre-commit locally, simply run

```bash
pre-commit install
```

in the root of the repository.
Pre-commit will automatically download all dependencies when it is run for the first time.

Alternatively, you can rely on the [pre-commit.ci][] service enabled on GitHub.
If you didn’t run `pre-commit` before pushing changes to GitHub it will automatically commit fixes to your pull request, or show an error message.

If pre-commit.ci added a commit on a branch you still have been working on locally, simply use

```bash
git pull --rebase
```

to integrate the changes into yours.
While the [pre-commit.ci][] is useful, we strongly encourage installing and running pre-commit locally first to understand its usage.

Finally, most editors have an _autoformat on save_ feature.
Consider enabling this option for [ruff][ruff-editors] and [biome][biome-editors].

[pre-commit]: https://pre-commit.com/
[pre-commit.ci]: https://pre-commit.ci/
[ruff-editors]: https://docs.astral.sh/ruff/integrations/
[biome-editors]: https://biomejs.dev/guides/integrate-in-editor/

(writing-tests)=

## Writing tests

This package uses [pytest][] for automated testing.
Please write {doc}`scanpy:dev/testing` for every function added to the package.

Most IDEs integrate with pytest and provide a GUI to run tests.
Just point yours to one of the environments returned by

```bash
hatch env create hatch-test  # create test environments for all supported versions
hatch env find hatch-test  # list all possible test environment paths
```

Alternatively, you can run all tests from the command line by executing

:::::{tabs}
::::{group-tab} Hatch

```bash
hatch test  # test with the highest supported Python version
# or
hatch test --all  # test with all supported Python versions
```

::::

::::{group-tab} uv

```bash
uv run pytest
```

::::

::::{group-tab} Pip

```bash
source .venv/bin/activate
pytest
```

::::
:::::

in the root of the repository.

[pytest]: https://docs.pytest.org/

### Continuous integration

Continuous integration via GitHub actions will automatically run the tests on all pull requests and test
against the minimum and maximum supported Python version.

Additionally, there’s a CI job that tests against pre-releases of all dependencies (if there are any).
The purpose of this check is to detect incompatibilities of new package versions early on and
gives you time to fix the issue or reach out to the developers of the dependency before the package
is released to a wider audience.

The CI job is defined in `.github/workflows/test.yaml`,
however the single point of truth for CI jobs is the Hatch test matrix defined in `pyproject.toml`.
This means that local testing via hatch and remote testing on CI tests against the same python versions and uses the same environments.

## Publishing a release

### Updating the version number

Before making a release, you need to update the version number in the `pyproject.toml` file.
Please adhere to [Semantic Versioning][semver], in brief

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> 1. MAJOR version when you make incompatible API changes,
> 2. MINOR version when you add functionality in a backwards compatible manner, and
> 3. PATCH version when you make backwards compatible bug fixes.
>
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Once you are done, commit and push your changes and navigate to the "Releases" page of this project on GitHub.
Specify `vX.X.X` as a tag name and create a release.
For more information, see [managing GitHub releases][].
This will automatically create a git tag and trigger a Github workflow that creates a release on [PyPI][].

[semver]: https://semver.org/
[managing GitHub releases]: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository
[pypi]: https://pypi.org/

## Writing documentation

Please write documentation for new or changed features and use-cases.
This project uses [sphinx][] with the following features:

- The [myst][] extension allows to write documentation in markdown/Markedly Structured Text
- [Numpy-style docstrings][numpydoc] (through the [napoloen][numpydoc-napoleon] extension).
- Jupyter notebooks as tutorials through [myst-nb][] (See [Tutorials with myst-nb](#tutorials-with-myst-nb-and-jupyter-notebooks))
- [sphinx-autodoc-typehints][], to automatically reference annotated input and output types
- Citations (like {cite:p}`Virshup_2023`) can be included with [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/)

See scanpy’s {doc}`scanpy:dev/documentation` for more information on how to write your own.

[sphinx]: https://www.sphinx-doc.org/en/master/
[myst]: https://myst-parser.readthedocs.io/en/latest/intro.html
[myst-nb]: https://myst-nb.readthedocs.io/en/latest/
[numpydoc-napoleon]: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[sphinx-autodoc-typehints]: https://github.com/tox-dev/sphinx-autodoc-typehints

### Tutorials with myst-nb and jupyter notebooks

The documentation is set-up to render jupyter notebooks stored in the `docs/notebooks` directory using [myst-nb][].
Currently, only notebooks in `.ipynb` format are supported that will be included with both their input and output cells.
It is your responsibility to update and re-run the notebook whenever necessary.

If you are interested in automatically running notebooks as part of the continuous integration,
please check out [this feature request][issue-render-notebooks] in the `cookiecutter-scverse` repository.

[issue-render-notebooks]: https://github.com/scverse/cookiecutter-scverse/issues/40

#### Hints

- If you refer to objects from other packages, please add an entry to `intersphinx_mapping` in `docs/conf.py`.
  Only if you do so can sphinx automatically create a link to the external documentation.
- If building the documentation fails because of a missing link that is outside your control,
  you can add an entry to the `nitpick_ignore` list in `docs/conf.py`

(docs-building)=

### Building the docs locally

:::::{tabs}
::::{group-tab} Hatch

```bash
hatch run docs:build
hatch run docs:open
```

::::

::::{group-tab} uv

```bash
cd docs
uv run sphinx-build -M html . _build -W
(xdg-)open _build/html/index.html
```

::::

::::{group-tab} Pip

```bash
source .venv/bin/activate
cd docs
sphinx-build -M html . _build -W
(xdg-)open _build/html/index.html
```

::::
:::::
