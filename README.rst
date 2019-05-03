|PyPI| |bioconda| |Docs| |Build Status|

.. |PyPI| image:: https://img.shields.io/pypi/v/scanpy.svg
   :target: https://pypi.org/project/scanpy
.. |bioconda| image:: https://img.shields.io/badge/bioconda-🐍-blue.svg
   :target: http://bioconda.github.io/recipes/scanpy/README.html
.. |Docs| image:: https://readthedocs.com/projects/icb-scanpy/badge/?version=latest
   :target: https://scanpy.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/theislab/scanpy.svg?branch=master
   :target: https://travis-ci.org/theislab/scanpy
..
   .. |Coverage| image:: https://codecov.io/gh/theislab/scanpy/branch/master/graph/badge.svg
      :target: https://codecov.io/gh/theislab/scanpy

Scanpy – Single-Cell Analysis in Python
=======================================

Scanpy is a scalable toolkit for analyzing single-cell gene expression data.
It includes preprocessing, visualization, clustering, pseudotime and trajectory
inference and differential expression testing. The Python-based implementation
efficiently deals with datasets of more than one million cells.

Read the documentation_.
If Scanpy is useful for your research, consider citing `Genome Biology (2018)`_.

.. _documentation: https://scanpy.readthedocs.io
.. _Genome Biology (2018): https://doi.org/10.1186/s13059-017-1382-0

I have added source code for weighted sampled data. I have already preprocessed data and found the top fea PC's and then input to scanpy
to found louvain coomunities , marker genes and later plot in different formats like dotPlot, violinPlot and heatmap. I have updated scanpy
for weighted sampled data where each row has its weight, but this support for clustering and plotting. We can further update sparse PCA as 
well to support weighted data points. 
