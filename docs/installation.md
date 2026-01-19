# Installation

To use `scanpy` from another project, install it using your favourite environment manager:

::::{tab-set}

:::{tab-item} Hatch (recommended)
:sync: hatch

Adding `scanpy[leiden]` to your dependencies is enough.
See below for how to use Scanpy’s {ref}`dev-install-instructions`.
:::

:::{tab-item} Pip
:sync: pip

If you prefer to exclusively use pip run:

```console
$ pip install 'scanpy[leiden]'
```
:::

:::{tab-item} Conda
:sync: conda

After installing installing e.g. [Miniconda][], run:

```console
$ conda install -c conda-forge scanpy python-igraph leidenalg
```

Pull Scanpy [from PyPI][]:

```console
$ pip install scanpy
```

[miniconda]: https://docs.anaconda.com/miniconda/miniconda-install/
[from pypi]: https://pypi.org/project/scanpy
:::

::::

If you use Hatch or pip, the extra `[leiden]` installs two packages that are needed for popular
parts of scanpy but aren't requirements: [igraph][] {cite:p}`Csardi2006` and [leiden][] {cite:p}`Traag2019`.
If you use conda, you should to add these dependencies to your environment individually.

[igraph]: https://python.igraph.org/en/stable/
[leiden]: https://leidenalg.readthedocs.io

(dev-install-instructions)=

## Development Version

To work with the latest version [on GitHub][]: clone the repository and `cd` into its root directory.

```console
$ gh repo clone scverse/scanpy
$ cd scanpy
```

::::{tab-set}

:::{tab-item} Hatch (recommended)
:sync: hatch

To use one of the predefined [Hatch environments][] in {file}`hatch.toml`,
run either `hatch test [args]` or `hatch run [env:]command [...args]`, e.g.:

```console
$ hatch test -p               # run tests in parallel
$ hatch run docs:build        # build docs
$ hatch run towncrier:create  # create changelog entry
```

[hatch environments]: https://hatch.pypa.io/latest/tutorials/environment/basic-usage/
:::

:::{tab-item} Pip
:sync: conda

If you are using `pip>=21.3`, an editable install can be made:

```console
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install --group=dev --group=test -e .
```
:::

:::{tab-item} Conda
:sync: conda

If you want to let `conda` handle the installations of dependencies, do:

```console
$ pipx install beni
$ beni pyproject.toml > environment.yml
$ conda env create -f environment.yml
$ conda activate scanpy
$ pip install --group=dev --group=test --group=doc -e .
```

For instructions on how to work with the code, see the {ref}`contribution guide <contribution-guide>`.
:::

::::

[on github]: https://github.com/scverse/scanpy

## Docker

If you're using [Docker][], you can use e.g. the image [gcfntnu/scanpy][] from Docker Hub.

[docker]: https://en.wikipedia.org/wiki/Docker_(software)
[gcfntnu/scanpy]: https://hub.docker.com/r/gcfntnu/scanpy

## Troubleshooting

If you get a `Permission denied` error, never use `sudo pip`. Instead, use virtual environments or:

```console
$ pip install --user scanpy
```
