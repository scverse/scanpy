**Plans for upcoming versions:**

- common file format for exchanging :class:`~scanpy.api.AnnData` with packages such as Seurat and SCDE
- better correction for confounders
- partial loading of data into memory

**Upcoming features for version 0.3.3**

- Export to `SPRING <https://github.com/AllonKleinLab/SPRING/>`_ [Weinreb17]_ for interactive visualization of data: `tutorial <https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export>`_, `docs <https://scanpy.readthedocs.io/en/latest/api/index.html>`_
  

**November 29, 2017**: version 0.3.2

- consistency updates, better logging, docs and bug fixes throughout
- Finding marker genes via :func:`~scanpy.api.pl.rank_genes_groups_violin` improved: `example <https://github.com/theislab/scanpy/issues/51>`_


**November 16, 2017**: version 0.3

- consistency updates, stability improvements
- :func:`~scanpy.api.AnnData` can be `concatenated <https://scanpy.readthedocs.io/en/latest/api/scanpy.api.AnnData.html>`_
- :func:`~scanpy.api.AnnData` is available as a `separate package <https://pypi.python.org/pypi/anndata/>`_
- results of approximate graph abstraction (AGA) are `simplified <https://github.com/theislab/graph_abstraction>`_


**October 25, 2017**: version 0.2.9

- initial release of `approximate graph abstraction (AGA) <https://github.com/theislab/graph_abstraction>`_.


**July 24, 2017**: version 0.2.1

Scanpy now includes preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing and simulation of gene regulatory networks. The implementation efficiently deals with datasets of more than one million cells.


**May 1, 2017**: version 0.1

Scanpy computationally outperforms the Cell Ranger R kit and allows reproducing most of Seurat's guided clustering tutorial.
