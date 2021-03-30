Installation
------------

Anaconda
~~~~~~~~
If you do not have a working installation of Python 3.6 (or later), consider
installing Miniconda_ (see `Installing Miniconda`_). Then run::

    conda install seaborn scikit-learn statsmodels numba pytables
    conda install -c conda-forge python-igraph leidenalg

Pull Scanpy from `PyPI <https://pypi.org/project/scanpy>`__ (consider using ``pip3`` to access Python 3)::

    pip install scanpy

.. _from PyPI: https://pypi.org/project/scanpy

PyPI only
~~~~~~~~~
If you prefer to exclusively use PyPI run::

    pip install 'scanpy[leiden]'

The extra `[leiden]` installs two packages that are needed for popular
parts of scanpy but aren't requirements: python-igraph_ [Csardi06]_ and leiden_ [Traag18]_.

.. _python-igraph: http://igraph.org/python/
.. _leiden: https://leidenalg.readthedocs.io

.. _dev-install-instructions:

Development Version
~~~~~~~~~~~~~~~~~~~
To work with the latest version `on GitHub`_: clone the repository and `cd` into
its root directory. To install using symbolic links (stay up to date with your
cloned version after you update with `git pull`) call::

    flit install -s --deps=develop  # from an activated venv or conda env
    # or
    flit install -s --deps=develop --python path/to/venv/bin/python

.. _on GitHub: https://github.com/theislab/scanpy

If you want to let conda_ handle the installations of dependencies, do::

    pip install beni
    beni pyproject.toml > environment.yml
    conda env create -f environment.yml
    conda activate scanpy
    flit install -s --deps=develop

On Windows, you might have to use `flit install --pth-file`
if you are not able to give yourself the `create symbolic links`_ privilege.

.. _create symbolic links: https://docs.microsoft.com/en-us/windows/security/threat-protection/security-policy-settings/create-symbolic-links

.. note::

    `pip install -e` still works, but may not in future versions.

Docker
~~~~~~
If you're using Docker_, you can use e.g. the image `gcfntnu/scanpy`_ from Docker Hub.

.. _Docker: https://en.wikipedia.org/wiki/Docker_(software)
.. _gcfntnu/scanpy: https://hub.docker.com/r/gcfntnu/scanpy
.. _bioconda: https://bioconda.github.io/

Troubleshooting
~~~~~~~~~~~~~~~
If you get a `Permission denied` error, never use `sudo pip`. Instead, use virtual environments or::

    pip install --user scanpy

**On MacOS**, if **not** using `conda`, you might need to install the C core of igraph via homebrew first

- `brew install igraph`
- If python-igraph still fails to install, see the question on `compiling igraph`_.
  Alternatively consider installing gcc via `brew install gcc --without-multilib`
  and exporting the required variables::

      export CC="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"
      export CXX="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"

  where `X` and `x` refers to the version of `gcc`;
  in my case, the path reads `/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6`.

**On Windows**, there also often problems installing compiled packages such as `igraph`,
but you can find precompiled packages on Christoph Gohlke’s `unofficial binaries`_.
Download those and install them using `pip install ./path/to/file.whl`

.. _compiling igraph: https://stackoverflow.com/q/29589696/247482
.. _unofficial binaries: https://www.lfd.uci.edu/~gohlke/pythonlibs/

.. _conda:

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~
After downloading Miniconda_, in a unix shell (Linux, Mac), run

.. code:: shell

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions.
Either reopen a new terminal or `source ~/.bashrc` on Linux/ `source ~/.bash_profile` on Mac.
The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html
