"""Run Diffusion maps using the adaptive anisotropic kernel
"""
from .. import logging as logg

def palantir( adata, **kargs ):
	
	"""
		Run Diffusion maps using the adaptive anisotropic kernel [Setty27]_.
		
		
		:param adata: :class:`~anndata.AnnData`, or Dataframe of cells X genes\n
		
		:param normalize: `bool` (default: `False`) provide argument
		                  for raw counts to normalize using palantir method, 
		                  `palantir.preprocess.normalize_counts`
		                  
		:param log_transform: `bool` (default: `False`) some datasets show 
		                  better signal in the log scale, applied using 
		                  `palantir.preprocess.log_transform`
		                  
		:return:
				  `data_df` - DataFrame of normalized adata\n
				  `pca_projections` - PCA projections of the data\n
				  `var_r` - PCA explained variance ratio\n
				  `dm_res` - Diffusion components, corresponding eigen values 
				  			 and the diffusion operator\n
				  `ms_data` - Multi scale data matrix\n
				  `tsne` - tSNE embedding of the data\n
				  `imp_df` - Imputed data matrix\n

	    Return objects will be pushed to adata.

	    Example
	    -------
	"""
	import numpy as np

	# logg.info('Palantir diffusion maps')

	class _wrapper_cls( object ):

		def __init__( self ,
								adata,
								func=None ,
								normalize = False,
								log_transform = False,
								copy = False
					):

			self.func = func
			self.adata = adata

			# load palantir
			self.__call__()
			# logg.info('palantir loaded ...')

			# load data and normalize if necessary
			self.preprocessing( self.adata,
								normalize = normalize,
								log_transform = log_transform )

			adata.obsm['palantir_norm_data'] = np.array(self.data_df)

			# Principal component analysis
			# logg.info('PCA in progress ...')
			self.pca_projections, self.var_r = self.palantir.utils.run_pca(self.data_df)
			adata.uns['palantir_pca'] = {}
			adata.uns['palantir_pca']['pca_projections'] = self.pca_projections
			adata.uns['palantir_pca']['variance_ratio'] = self.var_r

			# Diffusion maps
			# logg.info('Diffusion maps in progress ...')
			self.dm_res = self.palantir.utils.run_diffusion_maps(self.pca_projections)
			self.ms_data = self.palantir.utils.determine_multiscale_space(self.dm_res)
			adata.uns['palantir_diff_maps'] = {}
			adata.uns['palantir_diff_maps']['dm_res'] = self.dm_res
			adata.obsm['X_palantir_diffmap'] = np.array(self.ms_data)

			# tSNE visualization
			# logg.info('tSNE in progress ...')
			self.tsne = self.palantir.utils.run_tsne(self.ms_data)
			adata.obsm['X_palantir_tsne'] = np.array(self.tsne)

			# MAGIC imputation
			# logg.info('imputation in progress ...')
			self.imp_df = self.palantir.utils.run_magic_imputation(self.data_df, self.dm_res)
			adata.obsm['X_palantir_imputation'] = np.array(self.imp_df)

			# logg.info('End of processing, start plotting.')

		def __call__( self ):

			self.palantir = self.func()

		def preprocessing( self, data_df = None,
								 normalize = False,
								 log_transform = False):
			try:
				self.data_df = self.adata.to_df()
			except AttributeError:
				# assume the data is a cell X genes Dataframe
				pass
			print("data loaded ...")
			if normalize:
				try:
					self.data_df = self.palantir.preprocess.normalize_counts(self.data_df)
					print("data normalized ...")
				except AttributeError as e:
					raise AttributeError( "Missing Anndata: " + str(e) )
			if log_transform:
				try:
					self.data_df = self.palantir.preprocess.log_transform(self.data_df)
					print("data log transformed ...")
				except AttributeError as e:
					raise AttributeError( "Missing Anndata: " + str(e) )


	def wrapper_cls( adata, func=None , **kargs):
		if func:
			return _wrapper_cls( func )
		else:
			def wrapper( func ):
				return _wrapper_cls( adata, func , **kargs)
			return wrapper


	@wrapper_cls( adata , **kargs )
	def _run():
		import importlib
		try:
			palantir = importlib.import_module('palantir')
		except ImportError:
			raise ImportError(
				'\nplease install palantir: \n\n\t'
				'git clone git://github.com/dpeerlab/Palantir.git\n\t'
				'cd Palantir\n\t'
				'sudo -H pip3 install .\n')
		return palantir
	return _run
