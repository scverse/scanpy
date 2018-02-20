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

If you're using `Docker <https://en.wikipedia.org/wiki/Docker_(software)>`_, you
can setup your environment `here
<https://gist.github.com/flying-sheep/0e003ae3398dd543638955a55c031c8d>`_.

  
Trouble shooting
~~~~~~~~~~~~~~~~

If you do not have sudo rights (you get a ``Permission denied`` error)::

    pip install --user scanpy

**On MacOS**, you probably need to install the C core of igraph via homebrew first

- ``brew install igraph``
- If python-igraph still fails to install, see `here <https://stackoverflow.com/questions/29589696/problems-compiling-c-core-of-igraph-with-python-2-7-9-anaconda-2-2-0-on-mac-osx>`__ or consider installing gcc via ``brew install gcc --without-multilib`` and exporting ``export CC="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"; export CXX="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"``, where ``X`` and ``x`` refers to the version of ``gcc``; in my case, the path reads ``/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6``.

**On Windows**, you can follow these steps - thanks to Nikolaos-Kosmas Chlis! After installing Miniconda::
  
    conda install pandas cython
    conda install -c vtraag python-igraph
    conda install -c vtraag louvain
    pip install scanpy
      

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

After downloading Miniconda_, in a unix shell (Linux, Mac), run

.. code:: shell

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions. Either reopen a new terminal or ``source ~/.bashrc`` on Linux/ ``source ~/.bash_profile`` on Mac. The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html
