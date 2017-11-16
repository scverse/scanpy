**Plans for upcoming versions:**

- common file format for exchanging :class:`~scanpy.api.AnnData` with packages such as Seurat and SCDE
- better correction for confounders
- partial loading of data into memory


**November 16, 2017**: version 0.3

Consistency updates in the design of :class:`~scanpy.api.AnnData` and the code base. Everything is backwards compatible.

- `AnnData can now be concatenated <https://scanpy.readthedocs.io/en/stable/api/scanpy.api.AnnData.html>`_
- `AnnData is also available as a separate package <https://pypi.python.org/pypi/anndata/>`_
- `results of approximate graph abstraction (AGA) are simplified <https://github.com/theislab/graph_abstraction>`_
- stability fixes


**October 25, 2017**: version 0.2.9

Initial release of `approximate graph abstraction (AGA) <https://github.com/theislab/graph_abstraction>`_.


**July 24, 2017**: version 0.2.1

Scanpy now includes preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing and simulation of gene regulatory networks. The implementation efficiently deals with datasets of more than one million cells.


**May 1, 2017**: version 0.1

Scanpy computationally outperforms the Cell Ranger R kit and allows reproducing most of Seurat's guided clustering tutorial.
