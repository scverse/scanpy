def bbknn(adata, copy=False, **kwargs):
	'''
	Batch balanced KNN, altering the KNN procedure to identify each cell's top neighbours in
	each batch separately instead of the entire cell pool with no accounting for batch.
	Aligns batches in a quick and lightweight manner.
	For use in the scanpy workflow as an alternative to ``scanpi.api.pp.neighbors()``.
	
	For further details, check BBKNN's `GitHub <https://github.com/Teichlab/bbknn>`__ page 
	and [Park18]_. This docstring is accurate as of BBKNN version 1.2.0.
	
	Input
	-----
	adata : ``AnnData``
		Needs the PCA computed and stored in ``adata.obsm["X_pca"]``.
	batch_key : ``str``, optional (default: "batch")
		``adata.obs`` column name discriminating between your batches.
	neighbors_within_batch : ``int``, optional (default: 3)
		How many top neighbours to report for each batch; total number of neighbours 
		will be this number times the number of batches.
	n_pcs : ``int``, optional (default: 50)
		How many principal components to use in the analysis.
	trim : ``int`` or ``None``, optional (default: ``None``)
		If not ``None``, trim the neighbours of each cell to these many top connectivities.
		May help with population independence and improve the tidiness of clustering.
	scale_distance : ``bool``, optional (default: ``False``) 
		If ``True``, optionally lower the across-batch distances on a per-cell, per-batch basis to make
		the closest neighbour be closer to the furthest within-batch neighbour. 
		May help smooth out very severe batch effects with a risk of overly 
		connecting the cells. The exact algorithm is as follows:
		
		.. code-block:: python
		
			if min(corrected_batch) > max(original_batch):
				corrected_batch += max(original_batch) - min(corrected_batch) + np.std(corrected_batch)
	approx : ``bool``, optional (default: ``False``)
		If ``True``, use annoy's approximate neighbour finding. This results in a quicker run time 
		for large datasets at a risk of loss of independence of some of the populations. It should
		be noted that annoy's default metric of choice is "angular", which BBKNN overrides to
		"euclidean" from its own default metric setting.
	metric : ``str`` or ``sklearn.neighbors.DistanceMetric``, optional (default: "euclidean")
		What distance metric to use. If using ``approx=True``, the options are "euclidean",
		"angular", "manhattan" and "hamming". Otherwise, the options are "euclidean", 
		"manhattan", "chebyshev", or parameterised ``sklearn.neighbors.DistanceMetric`` 
		for "minkowski", "wminkowski", "seuclidean" or "mahalanobis".
		
		>>> from sklearn.neighbors import DistanceMetric
		>>> pass_this_as_metric = DistanceMetric.get_metric('minkowski',p=3)
	bandwidth : ``float``, optional (default: 1)
		``scanpy.neighbors.compute_connectivities_umap`` parameter, higher values result in a
		gentler slope of the connectivities exponentials (i.e. larger connectivity values being returned)
	local_connectivity : ``int``, optional (default: 1)
		``scanpy.neighbors.compute_connectivities_umap`` parameter, how many nearest neighbors of
		each cell are assumed to be fully connected (and given a connectivity value of 1)
	n_jobs : ``int`` or ``None``, optional (default: ``None``)
		Parallelise neighbour identification when using an Euclidean distance metric, 
		if ``None`` use all cores. Does nothing with a different metric.
	save_knn : ``bool``, optional (default: ``False``)
		If ``True``, save the indices of the nearest neighbours for each cell in ``adata.uns['bbknn']``.
	copy : ``bool``, optional (default: ``False``)
		If ``True``, return a copy instead of writing to the supplied adata.
	'''
	try:
		import bbknn
	except ImportError:
		raise ImportError('Install BBKNN: pip3 install bbknn')
	if copy:
		adata = bbknn.bbknn(adata, copy=copy, **kwargs)
		return adata
	else:
		bbknn.bbknn(adata, copy=copy, **kwargs)
