# Installation

## Anaconda

If you do not have a working installation of Python 3.6 (or later), consider
installing [Miniconda] (see [Installing Miniconda]). Then run:

```shell
conda install -c conda-forge scanpy python-igraph leidenalg
```

Pull Scanpy from [PyPI](https://pypi.org/project/scanpy) (consider using `pip3` to access Python 3):

```shell
pip install scanpy
```

## PyPI only

If you prefer to exclusively use PyPI run:

```shell
pip install 'scanpy[leiden]'
```

The extra `[leiden]` installs two packages that are needed for popular
parts of scanpy but aren't requirements: [igraph] {cite:p}`Csardi2006` and [leiden] {cite:p}`Traag2019`.

(dev-install-instructions)=

## Development Version

To work with the latest version [on GitHub]: clone the repository and `cd` into its root directory.

```shell
gh repo clone scverse/scanpy
cd scanpy
```

If you are using `pip>=21.3`, an editable install can be made:

```shell
pip install -e '.[dev,test]'
```

If you want to let [conda] handle the installations of dependencies, do:

```shell
pipx install beni
beni pyproject.toml > environment.yml
conda env create -f environment.yml
conda activate scanpy
pip install -e '.[dev,doc,test]'
```

For instructions on how to work with the code, see the {ref}`contribution guide <contribution-guide>`.

## Docker

If you're using [Docker], you can use e.g. the image [gcfntnu/scanpy] from Docker Hub.

## Troubleshooting

If you get a `Permission denied` error, never use `sudo pip`. Instead, use virtual environments or:

```shell
pip install --user scanpy
```

**On MacOS**, if **not** using `conda`, you might need to install the C core of igraph via homebrew first

- `brew install igraph`

- If igraph still fails to install, see the question on [compiling igraph].
  Alternatively consider installing gcc via `brew install gcc --without-multilib`
  and exporting the required variables:

  ```shell
  export CC="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"
  export CXX="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"
  ```

  where `X` and `x` refers to the version of `gcc`;
  in my case, the path reads `/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6`.

**On Windows**, there also often problems installing compiled packages such as `igraph`,
but you can find precompiled packages on Christoph Gohlke’s [unofficial binaries].
Download those and install them using `pip install ./path/to/file.whl`

(conda)=

## Installing Miniconda

After downloading [Miniconda], in a unix shell (Linux, Mac), run

```shell
cd DOWNLOAD_DIR
chmod +x Miniconda3-latest-VERSION.sh
./Miniconda3-latest-VERSION.sh
```

and accept all suggestions.
Either reopen a new terminal or `source ~/.bashrc` on Linux/ `source ~/.bash_profile` on Mac.
The whole process takes just a couple of minutes.

[bioconda]: https://bioconda.github.io/
[compiling igraph]: https://stackoverflow.com/q/29589696/247482
[create symbolic links]: https://docs.microsoft.com/en-us/windows/security/threat-protection/security-policy-settings/create-symbolic-links
[docker]: https://en.wikipedia.org/wiki/Docker_(software)
[from pypi]: https://pypi.org/project/scanpy
[gcfntnu/scanpy]: https://hub.docker.com/r/gcfntnu/scanpy
[leiden]: https://leidenalg.readthedocs.io
[miniconda]: https://docs.conda.io/projects/miniconda/en/latest/
[on github]: https://github.com/scverse/scanpy
[igraph]: https://python.igraph.org/en/stable/
[unofficial binaries]: https://www.lfd.uci.edu/~gohlke/pythonlibs/
