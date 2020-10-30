from anndata import AnnData
from typing import Collection, Tuple, Optional, Union
import numpy as np

from ... import logging as logg
from ... import preprocessing as pp

def scrub_doublets(
    adata: AnnData,
    sim_doublet_ratio: float = 2.0, 
    expected_doublet_rate: float = 0.05,
    stdev_doublet_rate: float =  0.02,
    synthetic_doublet_umi_subsampling: float = 1.0, 
    knn_dist_metric: str = 'euclidean', 
    normalize_variance: bool = True, 
    log_transform: bool = False,
    mean_center: bool = True,
    n_prin_comps: int = 30,
    use_approx_neighbors: bool = True, 
    get_doublet_neighbor_parents: bool = False, 
    n_neighbors: Optional[int] = None,
    threshold: Optional[float] = None,
    verbose: bool = True,
    copy: bool = False,
    random_state: int = 0
) -> Optional[Union[AnnData, Tuple[AnnData, AnnData], Tuple[AnnData, AnnData, AnnData]]]:
    """\
    Predict doublets using Scrublet [Wolock2019]_

    Predict cell doublets using a nearest-neighbor classifier of observed transcriptomes and simulated doublets. Works best if the input is a raw (unnormalized) counts matrix from a single sample or a collection of similar samples from the same experiment.
   
    This function is a wrapper around functions that pre-process using Scanpy and directly call functions of Scrublet(). You may also undertake your own preprocessing, simulate doublets with scanpy.external.pp.scrublet.simulate_doublets(), and run the core scrublet function scanpy.external.pp.scrublet.scrublet(). 

    .. note::
        More information and bug reports `here <https://github.com/swolock/scrublet>`__.
    
    Parameters
    ----------
    adata
        The annotated data matrix of shape ``n_obs`` × ``n_vars``. Rows
        correspond to cells and columns to genes. Expected to be un-normalised. 
    sim_doublet_ratio : float, optional (default: 2.0)
        Number of doublets to simulate relative to the number of observed 
        transcriptomes.
    expected_doublet_rate
        The estimated doublet rate for the experiment.
    stdev_doublet_rate
        Uncertainty in the expected doublet rate.
    synthetic_doublet_umi_subsampling
        Rate for sampling UMIs when creating synthetic doublets. If 1.0, each
        doublet is created by simply adding the UMI counts from two randomly
        sampled observed transcriptomes. For values less than 1, the UMI counts
        are added and then randomly sampled at the specified rate. 
    knn_dist_metric : str, optional (default: 'euclidean')
        Distance metric used when finding nearest neighbors. For list of
        valid values, see the documentation for annoy (if `use_approx_neighbors`
        is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
        is False).
    normalize_variance : bool, optional (default: True)
        If True, normalize the data such that each gene has a variance of 1.
        `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
        reduction, unless `mean_center` is True.
    log_transform
        Whether to use :func:``~scanpy.pp.log1p`` to log-transform the data prior to PCA.
    mean_center : bool, optional (default: True)
        If True, center the data such that each gene has a mean of 0.
        `sklearn.decomposition.PCA` will be used for dimensionality
        reduction.
    n_prin_comps : int, optional (default: 30)
        Number of principal components used to embed the transcriptomes prior
        to k-nearest-neighbor graph construction. 
    use_approx_neighbors : bool, optional (default: True)
        Use approximate nearest neighbor method (annoy) for the KNN 
        classifier.
    get_doublet_neighbor_parents : bool, optional (default: False)
        If True, return (in .uns) the parent transcriptomes that generated the
        doublet neighbors of each observed transcriptome. This information can
        be used to infer the cell states that generated a given doublet state. 
    n_neighbors
        Number of neighbors used to construct the KNN graph of observed
        transcriptomes and simulated doublets. If ``None``, this is
        automatically set to ``np.round(0.5 * np.sqrt(n_obs))``.
    threshold
        Doublet score threshold for calling a transcriptome a doublet. If
        `None`, this is set automatically by looking for the minimum between
        the two modes of the `doublet_scores_sim_` histogram. It is best
        practice to check the threshold visually using the
        `doublet_scores_sim_` histogram and/or based on co-localization of
        predicted doublets in a 2-D embedding.
    verbose : bool, optional (default: True)
    copy
        If ``True``, return a copy of the input ``adata`` with Scrublet results added. Otherwise, Scrublet results are added in place.
    random_state
        Initial state for doublet simulation and nearest neighbors.

    Returns
    -------
    adata : anndata.AnnData
        if ``copy=True`` it returns or else adds fields to ``adata``:
            ``.obs['doublet_score']``
                Doublet scores for each observed transcriptome
            
            ``.obs['predicted_doublets']``
                Boolean indicating predicted doublet status

            ``adata.uns['scrublet']['doublet_scores_sim']``
                Doublet scores for each simulated doublet transcriptome

            ``adata.uns['scrublet']['doublet_parents']`` 
                Pairs of ``.obs_names`` used to generate each simulated doublet transcriptome

            ``uns['scrublet']['parameters']``
                Dictionary of Scrublet parameters

    Examples
    --------

    """
    try:
        import scrublet as sl
    except ImportError:
        raise ImportError('Please install scrublet: `pip install scrublet` or `conda install scrublet`.')

    if copy:
        adata = adata.copy()

    start = logg.info('Running Scrublet')

    adata_obs = adata.copy()

    pp.filter_genes(adata_obs, min_cells = 3)
    pp.filter_cells(adata_obs, min_genes = 3)

    # Doublet simulation will be based on the un-normalised counts, but on the
    # selection of genes following normalisation and variability filtering. So
    # we need to save the raw and subset at the same time.

    adata_obs.layers['raw'] = adata_obs.X
    pp.normalize_total(adata_obs)

    # HVG process needs log'd data. If we're not using that downstream, then
    # copy logged data to new object and subset original object based on the
    # output.

    if log_transform:
        pp.log1p(adata_obs)
        pp.highly_variable_genes(adata_obs, subset = True)
    else:
        logged = pp.log1p(adata_obs, copy = True)
        hvg = pp.highly_variable_genes(logged)
        adata_obs = adata_obs[:, logged.var['highly_variable']]

    # Simulate the doublets based on the raw expressions from the normalised
    # and filtered object.

    adata_sim = simulate_doublets(
            adata_obs, 
            raw_layer = 'raw',
            sim_doublet_ratio = sim_doublet_ratio,
            synthetic_doublet_umi_subsampling = synthetic_doublet_umi_subsampling
    )

    # Now normalise simulated and observed in the same way

    pp.normalize_total(adata_obs, target_sum=1e6)
    pp.normalize_total(adata_sim, target_sum=1e6)

    adata_obs = scrublet(
            adata_obs = adata_obs,
            adata_sim = adata_sim,
            n_neighbors = n_neighbors,
            expected_doublet_rate = expected_doublet_rate,
            stdev_doublet_rate = stdev_doublet_rate,
            random_state = random_state
            )

    logg.info('    Scrublet finished', time = start)

    # Copy outcomes to input object from our processed version

    adata.obs['doublet_score'] = adata_obs.obs['doublet_score']
    adata.obs['predicted_doublet'] = adata_obs.obs['predicted_doublet']
    adata.uns['scrublet'] = adata_obs.uns['scrublet']

    if copy:
        return adata
    else:
        return None


def scrublet(
    adata_obs: AnnData,
    adata_sim: AnnData,
    n_neighbors: Optional[int] = None,
    expected_doublet_rate: float = 0.05,
    stdev_doublet_rate: float =  0.02,
    mean_center: bool = True,
    normalize_variance: bool = True, 
    n_prin_comps: int = 30,
    use_approx_neighbors: bool = True, 
    knn_dist_metric: str = 'euclidean', 
    get_doublet_neighbor_parents: bool = False, 
    random_state: int = 0, 
    verbose: bool = True,
) -> AnnData:
    """\
    Core function for predicting doublets using Scrublet [Wolock2019]_

    Predict cell doublets using a nearest-neighbor classifier of observed transcriptomes and simulated doublets. This is a wrapper around the core functions of `Scrublet <https://github.com/swolock/scrublet>`__ to allow for flexibility in applying Scanpy filtering operations upstream. Unless you know what you're doing you should use scrub_doublets().    
    
    .. note::
        More information and bug reports `here <https://github.com/swolock/scrublet>`__.
    
    Parameters
    ----------
    adata_obs
        The annotated data matrix of shape ``n_obs`` × ``n_vars``. Rows
        correspond to cells and columns to genes. Should be normalised with
        scanpy.pp.normalize_total() and filtered to include only highly
        variable genes.
    adata_sim
        Anndata object generated by
        sc.external.pp.scrublet.simulate_doublets(), with same number of vars
        as adata_obs. This should have been built from adata_obs after
        filtering genes and cells and selcting highly-variable genes..
    n_neighbors
        Number of neighbors used to construct the KNN graph of observed
        transcriptomes and simulated doublets. If ``None``, this is
        automatically set to ``np.round(0.5 * np.sqrt(n_obs))``.
    expected_doublet_rate
        The estimated doublet rate for the experiment.
    stdev_doublet_rate
        Uncertainty in the expected doublet rate.
    mean_center : bool, optional (default: True)
        If True, center the data such that each gene has a mean of 0.
        `sklearn.decomposition.PCA` will be used for dimensionality
        reduction.
    normalize_variance : bool, optional (default: True)
        If True, normalize the data such that each gene has a variance of 1.
        `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
        reduction, unless `mean_center` is True.
    n_prin_comps : int, optional (default: 30)
        Number of principal components used to embed the transcriptomes prior
        to k-nearest-neighbor graph construction. 
    use_approx_neighbors : bool, optional (default: True)
        Use approximate nearest neighbor method (annoy) for the KNN 
        classifier.
    knn_dist_metric : str, optional (default: 'euclidean')
        Distance metric used when finding nearest neighbors. For list of
        valid values, see the documentation for annoy (if `use_approx_neighbors`
        is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
        is False).
    get_doublet_neighbor_parents : bool, optional (default: False)
        If True, return the parent transcriptomes that generated the 
        doublet neighbors of each observed transcriptome. This information can 
        be used to infer the cell states that generated a given 
        doublet state. 
    random_state
        Initial state for doublet simulation and nearest neighbors.
    verbose : bool, optional (default: True)
        If True, print progress updates.

    Returns
    -------
    adata : anndata.AnnData
        if ``copy=True`` it returns or else adds fields to ``adata``:
            ``.obs['doublet_score']``
                Doublet scores for each observed transcriptome
            
            ``.obs['predicted_doublets']``
                Boolean indicating predicted doublet status

            ``adata.uns['scrublet']['doublet_scores_sim']``
                Doublet scores for each simulated doublet transcriptome

            ``adata.uns['scrublet']['doublet_parents']`` 
                Pairs of ``.obs_names`` used to generate each simulated doublet transcriptome

            ``uns['scrublet']['parameters']``
                Dictionary of Scrublet parameters

    Examples
    --------

    """
    try:
        import scrublet as sl
    except ImportError:
        raise ImportError('Please install scrublet: `pip install scrublet` or `conda install scrublet`.')

    # Estimate n_neighbors if not provided, and create scrublet object. 
    
    if n_neighbors is None:
        n_neighbors = int(round(0.5*np.sqrt(adata_obs.shape[0])))

    scrub = sl.Scrublet(
        adata_obs.X,
        n_neighbors = n_neighbors,
        expected_doublet_rate = expected_doublet_rate,
        stdev_doublet_rate = stdev_doublet_rate,
        random_state = random_state
        )

    scrub._E_obs_norm = adata_obs.X
    scrub._E_sim_norm = adata_sim.X
    scrub.doublet_parents_ = adata_sim.obsm['doublet_parents'] 
   
    # Call scrublet-specific preprocessing where specified

    if mean_center and normalize_variance:
        sl.pipeline_zscore(scrub)
    elif mean_center:
        sl.pipeline_mean_center(scrub)
    elif normalize_variance: 
        sl.pipeline_normalize_variance(scrub)

    # Do PCA. Scrublet fits to the observed matrix and decomposes both observed
    # and simulated based on that fit, so we'll just let it do its thing rather
    # than trying to use Scanpy's PCA wrapper of the same functions.
    
    if mean_center:
        logg.info('Embedding transcriptomes using PCA...')
        sl.pipeline_pca(scrub, n_prin_comps=n_prin_comps, random_state=scrub.random_state)
    else:
        logg.info('Embedding transcriptomes using Truncated SVD...')
        sl.pipeline_truncated_svd(scrub, n_prin_comps=n_prin_comps, random_state=scrub.random_state) 

    # Score the doublets

    scrub.calculate_doublet_scores(
        use_approx_neighbors=use_approx_neighbors,
        distance_metric=knn_dist_metric,
        get_doublet_neighbor_parents=get_doublet_neighbor_parents
    )

    # Actually call doublets

    scrub.call_doublets(verbose = verbose)

    # Store results in AnnData for return

    adata_obs.obs['doublet_score'] = scrub.doublet_scores_obs_
    adata_obs.obs['predicted_doublet'] = scrub.predicted_doublets_

    # Store doublet Scrublet metadata

    adata_obs.uns['scrublet'] = {}
    adata_obs.uns['scrublet']['threshold'] = scrub.threshold_
    adata_obs.uns['scrublet']['doublet_scores_sim'] = scrub.doublet_scores_sim_
    adata_obs.uns['scrublet']['doublet_parents'] = adata_sim.obsm['doublet_parents']
    adata_obs.uns['scrublet']['parameters'] = {
        'expected_doublet_rate': expected_doublet_rate,
        'sim_doublet_ratio': adata_sim.uns['scrublet']['parameters']['sim_doublet_ratio'],
        'n_neighbors': n_neighbors, 
        'random_state': random_state
    }

    if get_doublet_neighbor_parents:
        adata_obs.uns['scrublet']['doublet_neighbor_parents'] = scrub.doublet_neighbor_parents_
    
    return adata_obs

def simulate_doublets(
    adata: AnnData,
    raw_layer: str = 'raw',
    sim_doublet_ratio: float = 2.0,
    synthetic_doublet_umi_subsampling: float = 1.0,
    random_seed: int = 0
) -> AnnData:

    """\
    Simulate doublets by adding the counts of random observed transcriptome pairs.

    Arguments
    ---------
    adata
        The annotated data matrix of shape ``n_obs`` × ``n_vars``. Rows
        correspond to cells and columns to genes. Genes should have been
        filtered for expression and variability, and the object should contain
        raw expression of the same dimensions. 
        
    raw_layer : str, optional (default: 'raw')
        Layer of adata where raw values are stored, or 'X' if values are in .X. 
    
    sim_doublet_ratio : float, optional (default: None)
        Number of doublets to simulate relative to the number of observed 
        transcriptomes. If `None`, self.sim_doublet_ratio is used.

    synthetic_doublet_umi_subsampling : float, optional (defuault: 1.0) 
        Rate for sampling UMIs when creating synthetic doublets. If 1.0, 
        each doublet is created by simply adding the UMIs from two randomly 
        sampled observed transcriptomes. For values less than 1, the 
        UMI counts are added and then randomly sampled at the specified
        rate.

    Returns
    -------
    adata : anndata.AnnData with simulated doublets in .X
        if ``copy=True`` it returns or else adds fields to ``adata``:
            ``adata.uns['scrublet']['doublet_parents']`` 
                Pairs of ``.obs_names`` used to generate each simulated doublet transcriptome
            ``uns['scrublet']['parameters']``
                Dictionary of Scrublet parameters
    """
    try:
        import scrublet as sl
    except ImportError:
        raise ImportError('Please install scrublet: `pip install scrublet` or `conda install scrublet`.')

    if raw_layer == 'X':
        scrub = sl.Scrublet(adata.X)
    else:
        scrub = sl.Scrublet(adata.layers[raw_layer])
    
    scrub.simulate_doublets(
            sim_doublet_ratio = sim_doublet_ratio,
            synthetic_doublet_umi_subsampling =
            synthetic_doublet_umi_subsampling
    )
    
    adata_sim = AnnData(scrub._E_sim)
    adata_sim.obs['n_counts'] = scrub._total_counts_sim
    adata_sim.obsm['doublet_parents'] = scrub.doublet_parents_
    adata_sim.uns['scrublet'] = {}
    adata_sim.uns['scrublet']['parameters'] = {}
    adata_sim.uns['scrublet']['parameters']['sim_doublet_ratio'] = sim_doublet_ratio
    return adata_sim

def plot_histogram(adata, scale_hist_obs='log', scale_hist_sim='linear', fig_size = (8,3)):
    """\
    Plot histogram of doublet scores for observed transcriptomes and simulated doublets 

    The histogram for simulated doublets is useful for determining the correct doublet 
    score threshold. To set threshold to a new value, T, run call_doublets(threshold=T).
    """
    
    threshold = adata.uns['scrublet']['threshold']
    fig, axs = pl.subplots(1, 2, figsize=fig_size)

    ax = axs[0]
    ax.hist(adata.obs['doublet_score'], np.linspace(0, 1, 50), color='gray', linewidth=0, density=True)
    ax.set_yscale(scale_hist_obs)
    yl = ax.get_ylim()
    ax.set_ylim(yl)
    ax.plot(threshold * np.ones(2), yl, c='black', linewidth=1)
    ax.set_title('Observed transcriptomes')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Prob. density')

    ax = axs[1]
    ax.hist(adata.uns['scrublet']['doublet_scores_sim'], np.linspace(0, 1, 50), color='gray', linewidth=0, density=True)
    ax.set_yscale(scale_hist_sim)
    yl = ax.get_ylim()
    ax.set_ylim(yl)
    ax.plot(threshold * np.ones(2), yl, c='black', linewidth=1)
    ax.set_title('Simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Prob. density')

    fig.tight_layout()

    return fig, axs
