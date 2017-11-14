Installation
------------

If you do not have a working Python 3.5 or 3.6 installation, consider downloading and installing Miniconda_ (see `Installing Miniconda`_).

Get `releases on PyPI <https://pypi.python.org/pypi/scanpy>`__ via (consider using ``pip3`` to access Python 3)::

  pip install scanpy

To work with the latest version on `GitHub <https://github.com/theislab/scanpy>`__: clone the repository and ``cd`` into its root directory. To install using symbolic links (stay up to date with your cloned version after you update with ``git pull``) call::

    pip install -e .

Two further packages, needed for some of Scanpy's features, have not been automatically installed. Manually install them in this order

- `python-igraph <http://igraph.org/python/>`__ [Csardi06]_: ``pip install python-igraph`` 
- `louvain <https://github.com/vtraag/louvain-igraph>`__ [Traag17]_: ``pip install louvain>=0.6``

  
Trouble shooting
~~~~~~~~~~~~~~~~

If you do not have sudo rights (you get a ``Permission denied`` error)::

    pip install --user scanpy

**On MacOS**, you probably need to install the C core of igraph via homebrew first

- ``brew install igraph``
- If python-igraph still fails to install, see `here <https://stackoverflow.com/questions/29589696/problems-compiling-c-core-of-igraph-with-python-2-7-9-anaconda-2-2-0-on-mac-osx>`__ or consider installing gcc via ``brew install gcc --without-multilib`` and exporting ``export CC="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"; export CXX="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"``, where ``X`` and ``x`` refers to the version of ``gcc``; in my case, the path reads ``/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6``.

**On Windows**, you can follow these steps - thanks to Nikolaos-Kosmas Chlis!

- Step 0: Install Microsoft Visual C++ 10.0 standalone. Follow installation instructions `here <https://wiki.python.org/moin/WindowsCompilers#Microsoft_Visual_C.2B-.2B-_10.0_standalone:_Windows_SDK_7.1_.28x86.2C_x64.2C_ia64.29>`_. DO NOT ignore the first uninstallation step.

- Step 1: Download pycairo and igraph binaries from https://www.lfd.uci.edu/~gohlke/pythonlibs/::
    
    pycairo-1.15.2-cp35-cp35m-win_amd64.whl
    python_igraph-0.7.1.post6-cp35-none-win_amd64.whl
    
- Step 2: Create a conda environment called 'scanpyEnv' (or whatever you like) and install the rest via cmd::
  
    > conda create -n scanpyEnv python=3.5
    > activate scanpyEnv
    (scanpyEnv)> conda install pandas
    (scanpyEnv)> conda install cython
    (scanpyEnv)> pip install pycairo-1.15.2-cp35-cp35m-win_amd64.whl
    (scanpyEnv)> pip install python_igraph-0.7.1.post6-cp35-none-win_amd64.whl
    (scanpyEnv)> conda install -c vtraag louvain
    (scanpyEnv)> pip install scanpy
      

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

After downloading Miniconda_, in a unix shell (Linux, Mac), run

.. code:: shell

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions. Either reopen a new terminal or ``source ~/.bashrc`` on Linux/ ``source ~/.bash_profile`` on Mac. The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html
