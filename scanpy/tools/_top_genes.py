# Author: T. Callies
#
"""This modules provides all non-visualization tools for advanced gene ranking and exploration of genes
"""

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from .. import utils
from .. import logging as logg


def correlation_matrix(adata, name_list=None, groupby=None, group=None, n_genes=20, data='Complete', method='pearson', annotation_key=None):
    """Calculate correlation matrix.

        Calculate a correlation matrix for genes strored in sample annotation using :func:`~scanpy.api.tl.rank_genes_groups`.

        Parameters
        ----------
        adata : :class:`~anndata.AnnData`
            Annotated data matrix.
        name_list : list, optional (default: None)
            Takes a list of genes for which to calculate the correlation matrix
        groupby : `str`, optional (default: None)
            If no name list is passed, genes are selected from the
            results of rank_gene_groups. Then this is the key of the sample grouping to consider.
            Note that in this case also a group index has to be specified.
        group : `int`, optional (default: None)
            Group index for which the correlation matrix for top_ranked genes should be calculated.
            Currently only int is supported, will change very soon
        n_genes : `int`, optional (default: 20)
            For how many genes to calculate correlation matrix? If specified, cuts the name list
            (in whatever order it is passed).
        data : {'Complete', 'Group', 'Rest'}, optional (default: 'Complete')
            At the moment, this is only relevant for the case that name_list is drawn from rank_gene_groups results.
            If specified, collects mask for the called group and then takes only those cells specified.
            If 'Complete', calculate correlation using full data
            If 'Group', calculate correlation within the selected group.
            If 'Rest', calculate corrlation for everything except the group
        method : {‘pearson’, ‘kendall’, ‘spearman’} optional (default: 'pearson')
            Which kind of correlation coefficient to use
            pearson : standard correlation coefficient
            kendall : Kendall Tau correlation coefficient
            spearman : Spearman rank correlation
        annotation_key: String, optional (default: None)
            Allows to define the name of the anndata entry where results are stored.



    """

    # TODO: At the moment, only works for int identifiers

    ### If no genes are passed, selects ranked genes from sample annotation.
    ### At the moment, only calculate one table (Think about what comes next)
    if name_list is None:
        name_list = list()
        for j, k in enumerate(adata.uns['rank_genes_groups_gene_names']):
            if j >=n_genes:
                break
            name_list.append(adata.uns['rank_genes_groups_gene_names'][j][group])
    else:
        if len(name_list)>n_genes:
            name_list=name_list[0:n_genes]


    # If special method (later) , truncate
    adata_relevant = adata[:, name_list]
    # This line just makes group_mask access easier. Nothing else but 'all' will stand here.
    groups = 'all'
    if data == 'Complete' or groupby is None:
        if issparse(adata_relevant.X):
            Data_array = adata_relevant.X.todense()
        else:
            Data_array = adata_relevant.X
    else:
        # get group_mask
        groups_order, groups_masks = utils.select_groups(
            adata, groups, groupby)
        if data == 'Group':
            if issparse(adata_relevant.X):
                Data_array = adata_relevant.X[groups_masks[group], :].todense()
            else:
                Data_array = adata_relevant.X[groups_masks[group], :]
        elif data == 'Rest':
            if issparse(adata_relevant.X):
                Data_array = adata_relevant.X[~groups_masks[group], :].todense()
            else:
                Data_array = adata_relevant.X[~groups_masks[group], :]
        else:
            logg.error('data argument should be either <Complete> or <Group> or <Rest>' )

    # Distinguish between sparse and non-sparse data


    DF_array = pd.DataFrame(Data_array, columns=name_list)
    cor_table = DF_array.corr(method=method)
    if annotation_key is None:
        if groupby is None:
            adata.uns['Correlation_matrix'] = cor_table
        else:
            adata.uns['Correlation_matrix'+groupby+str(group)]=cor_table
    else:
        adata.uns[annotation_key] = cor_table




from sklearn import metrics
def ROC_AUC_analysis(adata,groupby,group=None, n_genes=100):
    """Calculate correlation matrix.

            Calculate a correlation matrix for genes strored in sample annotation using rank_genes_groups.py

            Parameters
            ----------
            adata : :class:`~anndata.AnnData`
                Annotated data matrix.
            groupby : `str`
                The key of the sample grouping to consider.
            group : `str`, int, optional (default: None)
                Group name or index for which the correlation matrix for top_ranked genes should be calculated.
                If no parameter is passed, ROC/AUC is calculated for all groups
            n_genes : `int`, optional (default: 100)
                For how many genes to calculate ROC and AUC. If no parameter is passed, calculation is done for
                all stored top ranked genes.

        """
    if group is None:
        pass
        # TODO: Loop over all groups instead of just taking one.

    # Assume group takes an int value for one group for the moment.
    name_list = list()
    for j, k in enumerate(adata.uns['rank_genes_groups_gene_names']):
        if j >= n_genes:
            break
        name_list.append(adata.uns['rank_genes_groups_gene_names'][j][group])


    # TODO: For the moment, see that everything works for comparison against the rest. Resolve issues later.
    groups = 'all'
    groups_order, groups_masks = utils.select_groups(
        adata, groups, groupby)

    # Use usual convention, better for looping later.
    imask = group
    mask = groups_masks[group]

    # TODO: Allow for sample weighting requires better mask access... later

    # We store calculated data in dict, access it via dict to dict. Check if this is the best way.
    fpr={}
    tpr={}
    thresholds={}
    roc_auc={}
    y_true=mask
    for i, j in enumerate(name_list):
        vec=adata[:,[j]].X
        if issparse(vec):
            y_score = vec.todense()
        else:
            y_score = vec

        fpr[name_list[i]], tpr[name_list[i]], thresholds[name_list[i]] = metrics.roc_curve(y_true, y_score, pos_label=None, sample_weight=None, drop_intermediate=False)
        roc_auc[name_list[i]]=metrics.auc(fpr[name_list[i]],tpr[name_list[i]])
    adata.uns['ROCfpr' +groupby+ str(group)] = fpr
    adata.uns['ROCtpr' +groupby+ str(group)] = tpr
    adata.uns['ROCthresholds' +groupby+ str(group)] = thresholds
    adata.uns['ROC_AUC' + groupby + str(group)] = roc_auc

def subsampled_estimates(mask, mask_rest=None, precision=0.01, probability=0.99):
    ## Simple method that can be called by rank_gene_group. It uses masks that have been passed to the function and
    ## calculates how much has to be subsampled in order to reach a certain precision with a certain probability
    ## Then it subsamples for mask, mask rest
    ## Since convergence speed varies, we take the slower one, i.e. the variance. This might have future speed-up
    ## potential
    if mask_rest is None:
        mask_rest=~mask
    # TODO: DO precision calculation for mean variance shared


    # TODO: Subsample

def dominated_ROC_elimination(adata,grouby):
    ## This tool has the purpose to take a set of genes (possibly already pre-selected) and analyze AUC.
    ## Those and only those are eliminated who are dominated completely
    ## TODO: Potentially (But not till tomorrow), this can be adapted to only consider the AUC in the given
    ## TODO: optimization frame
    pass


def _gene_preselection(adata,mask,thresholds):
    ## This tool serves to
    ## It is not thought to be addressed directly but rather using rank_genes_group or ROC analysis or comparable
    ## TODO: Pass back a truncated adata object with only those genes that fullfill thresholding criterias
    ## This function should be accessible by both rank_genes_groups and ROC_curve analysis
    pass
