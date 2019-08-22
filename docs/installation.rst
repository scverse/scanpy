Installation
------------

Anaconda
~~~~~~~~

If you do not have a working Python 3.5 or 3.6 installation, consider
installing Miniconda_ (see `Installing Miniconda`_). Then run::

    conda install seaborn scikit-learn statsmodels numba pytables
    conda install -c conda-forge python-igraph louvain

Pull Scanpy from `PyPI <https://pypi.org/project/scanpy>`__ (consider
using ``pip3`` to access Python 3)::

    pip install scanpy

PyPI only
~~~~~~~~~

If you prefer to exclusively use PyPI run::

    pip install scanpy[louvain]

or::
    
    pip install scanpy python-igraph louvain

The extra ``[louvain]`` installs two packages that are needed for
parts of scanpy but aren't requirements: `python-igraph
<http://igraph.org/python/>`__ [Csardi06]_ and `louvain
<https://github.com/vtraag/louvain-igraph>`__ [Traag17]_.

Bioconda
~~~~~~~~

Using bioconda_, simply run::

    conda install -c bioconda scanpy

Development Version
~~~~~~~~~~~~~~~~~~~

To work with the latest version on `GitHub
<https://github.com/theislab/scanpy>`__: clone the repository and ``cd`` into
its root directory. To install using symbolic links (stay up to date with your
cloned version after you update with ``git pull``) call::

    pip install -e .

Docker
~~~~~~

If you're using Docker_, you can use the minimal `fastgenomics/scanpy`_ image from the Docker Hub.

.. _Docker: https://en.wikipedia.org/wiki/Docker_(software)
.. _fastgenomics/scanpy: https://hub.docker.com/r/fastgenomics/scanpy
.. _bioconda: https://bioconda.github.io/

Troubleshooting
~~~~~~~~~~~~~~~

If you get a `Permission denied` error, never use `sudo pip`. Instead, use virtual environments or::

    pip install --user scanpy

**On MacOS**, if **not** using `conda`, you might need to install the C core of igraph via homebrew first

- ``brew install igraph``
- If python-igraph still fails to install, see `here <https://stackoverflow.com/questions/29589696/problems-compiling-c-core-of-igraph-with-python-2-7-9-anaconda-2-2-0-on-mac-osx>`__ or consider installing gcc via ``brew install gcc --without-multilib`` and exporting ``export CC="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"; export CXX="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"``, where ``X`` and ``x`` refers to the version of ``gcc``; in my case, the path reads ``/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6``.

**On Windows**, there also often problems installing compiled packages such as `igraph`, but you can find precompiled packages on <https://www.lfd.uci.edu/~gohlke/pythonlibs/>. Download those and install them using `pip install ./path/to/file.whl`

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

After downloading Miniconda_, in a unix shell (Linux, Mac), run

.. code:: shell

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions. Either reopen a new terminal or ``source ~/.bashrc`` on Linux/ ``source ~/.bash_profile`` on Mac. The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html
