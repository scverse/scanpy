from typing import Iterable

from anndata import AnnData

from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
import numpy as np
import itertools
import pandas as pd


def _demultiplex_per_barcode(x):
    """
    Approximates method in https://www.nature.com/articles/s41592-019-0433-8
    MULTI-seq: sample multiplexing for single-cell RNA
    sequencing using lipid-tagged indices
    """
    percentiles = np.percentile(x, q=np.arange(.01, 99.9, step=.1))
    kde = gaussian_kde(x)
    dist = kde(percentiles)
    maxima = argrelextrema(dist, np.greater_equal)[0]
    min_maxima = np.min(maxima)
    max_maxima = np.max(maxima)
    dist[min_maxima], dist[max_maxima]
    intermediate_thresh = (percentiles[min_maxima] + percentiles[max_maxima]) / 2
    all_low_maxima = maxima[np.where(percentiles[maxima] < intermediate_thresh)[0]]
    low_maxima = all_low_maxima[np.argmax(dist[all_low_maxima])]
    info = []
    for q in range(1, 100):
        threshold = np.percentile([percentiles[low_maxima],
                                   percentiles[max_maxima]], q=q)
        info.append([x >= threshold, threshold])
    return info, [percentiles, dist]


def demultiplex_hashing(adata: AnnData,
                        percentile_candidates: Iterable[int] = None,
                        inplace: bool = True):
    """
    Demultiplexing cell hashing data

    Params
    ------
    adata
        The annotated data matrix of shape ``n_obs`` × ``n_vars``. Rows correspond
        to cells and columns to genes.
    percentile_candidates
        Choose exact percentiles to evaluate thresholding at

    inplace
        Whether to update ``adata`` or return dictionary with normalized copies of
        ``adata.X`` and ``adata.layers``.

    Returns
    -------
    Returns `adata.obs` with columns "ID" and "CLASSIFICATION" which identify
    which cell hashing barcodes a cell has
    """
    barcodes_hits = []
    number_of_barcodes = adata.X.shape[1]
    dists_info = []
    for bc_idx in range(number_of_barcodes):
        x = adata[:, bc_idx].X
        info, dist_info = _demultiplex_per_barcode(x)
        dists_info.append(dist_info)
        barcodes_hits.append(info)
    number_of_singlets = []
    thresh_results = np.array(barcodes_hits)[:, :, 0]
    if percentile_candidates is not None:
        all_products = list(itertools.product(percentile_candidates,
                                              repeat=number_of_barcodes))
    else:
        all_products = [[q] * number_of_barcodes for q in range(1, 99)]
    for x in all_products:
        subsets = thresh_results[np.arange(number_of_barcodes), x]
        subsets = np.array([np.array(subset) for subset in subsets])
        number_of_singlets.append(np.sum(np.sum(subsets.astype(int),
                                                axis=0) == 1))
    percentiles_to_use = all_products[np.argmax(number_of_singlets)]
    thresholds_to_use = np.array(barcodes_hits)[np.arange(number_of_barcodes),
                                                percentiles_to_use, 1]

    barcode_hits = thresholds_to_use < adata.X[:]
    barcode_hit_indices = np.where(barcode_hits)
    cell_index_to_barcode_sub_ind = {idx:
                                     np.where(barcode_hit_indices[0] == hit)[0]
                                     for idx, hit in
                                     enumerate(range(adata.shape[0]))}

    cell_index_catergory_group = []
    for k, v in cell_index_to_barcode_sub_ind.items():
            if v is None:
                continue
            barcode_groups = \
                np.ravel(adata.var_names[barcode_hit_indices[1][v]])
            if len(barcode_groups) == 0:
                category = "Negative"
                group = "Negative"
            elif len(barcode_groups) == 1:
                category = barcode_groups[0]
                group = barcode_groups[0]
            else:
                category = "Doublet"
                group = "_".join(barcode_groups)
            cell_index_catergory_group.append([k, category, group])
    cell_index_catergory_group_df = pd.DataFrame(cell_index_catergory_group)[[1,2]]
    adata.obs["ID"] = cell_index_catergory_group_df[1].astype(str).tolist()
    adata.obs["CLASSIFICATION"] = cell_index_catergory_group_df[2].astype(str).tolist()
    return adata if not inplace else None

