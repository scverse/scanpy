Ecosystem
=========

.. role:: small
.. role:: smaller

With *ecosystem*, we mean single-cell related tools that operate on :class:`~anndata.AnnData`. Here, we list some that do not have an interface in the :doc:`external API <external/index>`.


Viewers
-------

Interactive manifold viewers.

* `cellxgene <https://github.com/chanzuckerberg/cellxgene>`__ via direct reading of `.h5ad` :small:`CZI`
* `cirrocumulus <https://cirrocumulus.readthedocs.io/>`__ via direct reading of `.h5ad` :small:`Broad Inst.`
* `cell browser <https://cells.ucsc.edu/>`__ via exporing through :func:`~scanpy.external.exporting.cellbrowser` :small:`UCSC`
* `SPRING <https://github.com/AllonKleinLab/SPRING>`__ via exporting through :func:`~scanpy.external.exporting.spring_project` :small:`Harvard Med`


Portals
-------

* the `Gene Expression Analysis Resource <https://umgear.org/>`__ :small:`U Maryland`
* the `Galaxy Project <https://humancellatlas.usegalaxy.eu>`__ for the Human Cell Atlas `[tweet] <https://twitter.com/ExpressionAtlas/status/1151797848469626881>`__ :small:`U Freiburg`
* the `Expression Atlas <https://www.ebi.ac.uk/gxa/sc/help.html>`__ :small:`EMBL-EBI`


RNA velocity
------------

* `scVelo <https://scvelo.org>`__ :small:`Helmholtz Munich`


Fate mapping
------------

* `CellRank <http://cellrank.org>`__ :small:`Helmholtz Munich`

    | CellRank is a toolkit to uncover cellular dynamics based on scRNA-seq data with
      RNA velocity annotation by detecting initial and terminal populations, inferring
      fate potentials and uncovering gene expression trends towards specific
      terminal populations.


Differential expression
-----------------------

* `diffxpy <https://github.com/theislab/diffxpy>`__ :small:`Helmholtz Munich`


Data integration
----------------

* `scanaroma <https://github.com/brianhie/scanorama>`__ :small:`MIT`


Modeling perturbations
----------------------

* `scGen <https://github.com/theislab/scgen>`__ / `trVAE <https://github.com/theislab/trvae>`__ :small:`Helmholtz Munich`


scvi-tools
----------

* `scvi-tools <https://github.com/YosefLab/scvi-tools>`__ :small:`Berkeley`

    | scvi-tools hosts deep generative models (DGM) for end-to-end analysis of single-cell
      omics data (e.g., scVI, scANVI, totalVI). It also contains several primitives to build novel DGMs.
      

Adaptive immune receptor repertoire (AIRR)
------------------------------------------

* `scirpy <https://github.com/icbi-lab/scirpy>`__ :small:`Medical University of Innsbruck`
    
    | scirpy is a scanpy extension to expore single-cell T-cell receptor (TCR) and B-cell receptor (BCR) repertoires. 
