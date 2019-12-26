Ecosystem
=========

Tools that operate on :class:`~anndata.AnnData` but do not have an interface in the :doc:`external API <external/index>`.

Viewers
-------

Interactive manifold viewers that directly read anndata's `.h5ad` or have an export in :doc:`external API <external/index>`

* `cellxgene <https://github.com/chanzuckerberg/cellxgene>`__ [CZI] via direct reading of `.h5ad`
* `cell browser <https://cells.ucsc.edu/>`__ [UCSC] via :func:`~scanpy.external.exporting.cellbrowser`
* `SPRING <https://github.com/AllonKleinLab/SPRING>`__ [Allon Klein lab] via :func:`~scanpy.external.exporting.spring_project`

Portals
-------

* The `Gene Expression Analysis Resource <https://umgear.org/>`__ [U Maryland]
* The `Galaxy Project <https://humancellatlas.usegalaxy.eu>`__ for Human Cell Atlas data [U Freiburg] `[tweet] <https://twitter.com/ExpressionAtlas/status/1151797848469626881>`__
* The EMBL-EBI `Expression Atlas <https://www.ebi.ac.uk/gxa/sc/help.html>`__


RNA velocity
------------

* `scvelo <https://scvelo.org>`__ [Helmholtz Munich]

Differential expression
-----------------------

* `diffxpy <https://github.com/theislab/diffxpy>`__ [Helmholtz Munich]

Integrating datasets and modeling data
--------------------------------------

* `scVI <https://github.com/YosefLab/scVI>`__ [Yosef lab]
* `scanaroma <https://github.com/brianhie/scanorama>`__ [Berger lab]

Modeling perturbations
----------------------

* `scgen <https://github.com/theislab/scgen>`__ (`trvae <https://github.com/theislab/trvae>`__) [Helmholtz Munich]
