# Author: T. Callies
#
"""This modules provides all visualization tools for advanced gene ranking and exploration of genes. They
are captured here and accessed through the standard function call sc.pl.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import issparse
from scipy.stats import rankdata
from .. import utils
from .. import logging as logg
from ..preprocessing import simple
import matplotlib.cm as cm
import pandas as pd

def correlation_matrix(adata,groupby=None ,group=None, corr_matrix=None, annotation_key=None):
    """Plot correlation matrix.

            Plot a correlation matrix for genes strored in sample annotation using rank_genes_groups.py

            Parameters
            ----------
            adata : :class:`~anndata.AnnData`
                Annotated data matrix.
            groupby : `str`, optional (default: None)
                If specified, searches data_annotation for correlation_matrix+groupby+str(group)
            group : int
                Identifier of the group (necessary if and only if groupby is also specified)
            corr_matrix : DataFrame, optional (default: None)
                Correlation matrix as a DataFrame (annotated axis) that can be transferred manually if wanted
            annotation_key: `str`, optional (default: None)
                If specified, looks in data annotation for this key.

        """

    # TODO: At the moment, noly works for int identifiers

    if corr_matrix is None:
        # This will produce an error if he annotation doesn't exist, which is okay
        if annotation_key is None:
            if groupby is None:
                corr_matrix = adata.uns['Correlation_matrix']
            else:
                corr_matrix= adata.uns['Correlation_matrix' + groupby+ str(group)]
            # Throws error if does not exist
        else:
            corr_matrix = adata.uns[annotation_key]

    # Set up mask
    mask = np.zeros_like(corr_matrix, dtype=np.bool)
    di = np.diag_indices(len(corr_matrix.axes[0]))
    mask[di] = True

    f, ax = plt.subplots(figsize=(11, 9))

    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    sns.heatmap(corr_matrix, mask=mask, cmap=cmap,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})
    if annotation_key is None:
        if groupby is None:
            plt.title('Correlation Matrix')
        else:
            plt.title('Correlation Matrix for Group ' + str(group) + "in " + groupby)
    else:
        plt.title('Correlation Matrix for' + annotation_key)
    plt.show()

def exploratory_rank_analysis(adata, groupby, x='inflation', y='mean', groups='all', n=100, special_markers=None,
                              coloring='scores',annotate=False):
    """Plot scatterplots for various gene_characteristics.

                This is a visualization tools that helps to find significant markers and get a better understanding of
                the underlying data.

                Parameters
                ----------
                adata : :class:`~anndata.AnnData`
                    Annotated data matrix.
                groupby : `str`
                    The key of the sample grouping to consider.
                x : 'str'
                    x-axis labelling for plots
                y : 'str'
                    y-axis labelling for plots
                groups : `str`, `list`, optional (default: `'all'`)
                    Subset of groups, e.g. `['g1', 'g2', 'g3']`, to which comparison shall
                    be restricted. If not passed, a ranking will be generated for all
                    groups.
                n : `int`, optional (default: 100)
                    Number of datapoints in the scatterplot. If less are available, use all that are available
                special_markers: 'dict', optional (default: None)
                    If provided, this should be a dict containing a list of gene names for each group in groupby.
                    Special marked genes are highlighted in the visualization
                coloring : {'scores', 'absolute'}, optional (default: 'scores')
                    Rank either according to Scores, or to absolute test-statistic value.
                    In either case, results are scaled so as to guarantee sufficient contrast.
                annotate: bool, optional (default: False)
                    If set to TRUE, annotate each scatterpoint with its name. Advisable only for small
                    number of data points.
            """

    # TODO: Check closely what of the below actions can be skipped and whats necessary


    n_groups = 0
    for i, j in enumerate(adata.uns['rank_genes_groups_gene_names'][0]):
        n_groups = n_groups + 1
    # Get group masks
    # TODO: Generalize. At the moment, only groups='all' works
    groups_order, groups_masks = utils.select_groups(
        adata, groups, groupby)
    # Create figure:
    n_rows = int(n_groups / 4) + 1
    n_cols = 4
    # For each group, get right genes (can be handled by intern function?)
    plt.figure(figsize=(24, 16))
    for imask, mask in enumerate(groups_masks):
        score_list = list()
        name_list = list()
        special_markers_indices = list()
        # Note: No duplicates in each group
        for j, k in enumerate(adata.uns['rank_genes_groups_gene_scores']):
            # Make sure only first n datapoints are used
            if j >= n:
                break
            score_list.append(k[imask])
            name_list.append(adata.uns['rank_genes_groups_gene_names'][j][imask])
            # Inefficient and not generalizable but works: Check added index if in list of specially_marked_genes
            # TODO: Speed up if becomes a time issue
            if special_markers is None:
                pass
            elif adata.uns['rank_genes_groups_gene_names'][j][imask] in special_markers[imask]:
                special_markers_indices.append(len(name_list) - 1)
            else:
                pass

        ### Get all the key figures
        # make things faster by calculating only what is required for plot
        mask_rest = ~mask

        # Get rate of expression
        rate_group = _zero_inflation_estimate(adata[:, name_list], mask)
        rate_rest = _zero_inflation_estimate(adata[:, name_list], mask_rest)
        if (x in {'full_mean_group', 'tail_mean_group','full_mean_difference',
                 'tail_mean_difference'} or y in {'full_mean_group', 'tail_mean_group','full_mean_difference',
                 'tail_mean_difference'}):
            means_group = _tail_mean_estimate(adata[:, name_list], mask)
        if (x in {'full_mean_rest', 'tail_mean_rest', 'full_mean_difference',
                  'tail_mean_difference'} or y in {'full_mean_rest', 'tail_mean_rest', 'full_mean_difference',
                                                   'tail_mean_difference'}):
            means_rest = _tail_mean_estimate(adata[:, name_list], mask_rest)

        if (x == 'tail_var_group' or y == 'tail_var_group'):
            # Get tail variance of expression
            var_group = _tail_var_estimate(adata[:, name_list], mask)
        if (x == 'tail_var_rest' or y == 'tail_var_rest'):
            var_rest = _tail_var_estimate(adata[:, name_list], mask_rest)
        if (x == 'CDR' or y == 'CDR'):
            # Get CDR: Need to give full adata object, since we need to count everything
            CDR = _Avg_CDR(adata, mask, name_list, model='rough', n_genes=None)
        if (x == 'full_var_group' or y == 'full_var_group'):
            # Slice first appropriately:
            adata_relevant = adata[:, name_list]
            exp, full_var_group = simple._get_mean_var(adata_relevant.X[mask])
        if (x == 'full_var_rest' or y == 'full_var_rest'):
            # Slice first appropriately:
            adata_relevant = adata[:, name_list]
            exp_rest, full_var_rest = simple._get_mean_var(adata_relevant.X[mask_rest])

        ### Prepare for coloring
        # get colored scatterplot
        # For coloring, get max score value, normalize (0,1)
        # Depending on whether normalization should be scale-invariant or only rank-invariant, do the following
        if coloring == 'scores':
            score_list = score_list / max(score_list)
            colors = cm.jet(score_list)
        elif coloring == 'absolute':
            color_list = rankdata(score_list)
            max_values = max(color_list)
            colors = cm.jet(color_list / max_values)
            # Identify true markers distinctly by using different size.
        else:
            logg.error('coloring should be either <socres> or <absolute>')
        s = 20 * np.ones(len(score_list))
        # This works for numpy access (not for normal lists though)
        s[special_markers_indices] = 100
        # In future, build method to mark top genes specially

        ### Actually do the plotting: Looping is inefficient and lengthy, but clear style
        # Potential values for x, y: 'mean' ('full' or 'tail'), 'tail_variance', 'inflation', 'CDR',
        # tail_variance_rest, Score (Just the ranking as given by test-statistic), 'full_var', 'full_var_rest'

        if x == 'expression_rate_difference':
            x_plot = rate_group - rate_rest
        elif x == 'expression_rate_group':
            x_plot = rate_group
        elif x == 'expression_rate_rest':
            x_plot = rate_rest
        elif x == 'Score':
            x_plot = score_list
        elif x == 'full_mean_difference':
            x_plot = means_group*rate_group-means_rest*rate_rest
        elif x == 'full_mean_group':
            x_plot = means_group*rate_group
        elif x == 'full_mean_rest':
            x_plot = means_rest*rate_rest
        elif x == 'tail_mean_difference':
            x_plot = means_group-means_rest
        elif x == 'tail_mean_group':
            x_plot = means_group
        elif x == 'tail_mean_rest':
            x_plot = means_rest
        elif x == 'tail_var_group':
            x_plot = var_group
        elif x == 'tail_var_rest':
            x_plot = var_rest
        elif x == 'full_var_group':
            x_plot = full_var_group
        elif x == 'full_var_rest':
            x_plot = full_var_rest
        elif x == 'CDR':
            x_plot = CDR
        else:
            logg.error('No accepted input. Check function documentation to get an overview over all inputs')

        if y == 'expression_rate_difference':
            y_plot = rate_group - rate_rest
        elif y == 'expression_rate_group':
            y_plot = rate_group
        elif y == 'expression_rate_rest':
            y_plot = rate_rest
        elif y == 'Score':
            y_plot = score_list
        elif y == 'full_mean_difference':
            y_plot = means_group*rate_group-means_rest*rate_rest
        elif y == 'full_mean_group':
            y_plot = means_group*rate_group
        elif y == 'full_mean_rest':
            y_plot = means_rest*rate_rest
        elif y == 'tail_mean_difference':
            y_plot = means_group-means_rest
        elif y == 'tail_mean_group':
            y_plot = means_group
        elif y == 'tail_mean_rest':
            y_plot = means_rest
        elif y == 'tail_var_group':
            y_plot = var_group
        elif y == 'tail_var_rest':
            y_plot = var_rest
        elif y == 'full_var_group':
            y_plot = full_var_group
        elif y == 'full_var_rest':
            y_plot = full_var_rest
        elif y == 'CDR':
            y_plot = CDR
        else:
            logg.error('No accepted input. Check function documentation to get an overview over all inputs')

        plt.subplot(n_rows, n_cols, imask + 1)
        plt.xlabel(x)
        plt.ylabel(y)

        # To make different scalings easier to compare, we set fixed limits for the case that x,y are e
        # expression rates
        if (x in {'expression_rate_difference', 'expression_rate_group', 'expression_rate_rest'} and
                    y in {'expression_rate_difference', 'expression_rate_group', 'expression_rate_rest'}):
            plt.xlim(0, 1)
            plt.ylim(0, 1)
        plt.scatter(x_plot, y_plot, color=colors, s=s)
        if annotate is True:
            for i, txt in enumerate(name_list):
                plt.annotate(txt, (x_plot[i], y_plot[i]))
    plt.show()

def _zero_inflation_estimate(adata, mask, model='rough'):
    # Method ZINB will be implemented soon
    if model not in {'rough', 'zinb'}:
        model = 'rough'
        logg.warn('Model should be either rough or zinb (zero-inflated negative binomial)')
    if adata.X.shape is int:
        X = adata.X[mask]
    else:
        X = adata.X[mask, :]
    n_cells = X.shape[0]
    if model == 'rough':
        if issparse(X):
            return X.getnnz(axis=0) / n_cells
        else:
            return np.count_nonzero(X, axis=0) / n_cells
    else:
        # Method for ZINB will be included soon
        if issparse(X):
            return 0
        else:
            return 0


def _tail_mean_estimate(adata, mask, model='rough'):
    # Method ZINB will be implemented soon
    if model not in {'rough', 'zinb'}:
        model = 'rough'
        logg.warn('Model should be either rough or zinb (zero-inflated negative binomial)')
    X = adata.X[mask, :]
    n_cells = X.shape[0]
    n_genes = X.shape[1]
    means = np.zeros((n_genes,))
    if model == 'rough':
        if issparse(X):
            n_nonzero_elements = X.getnnz(axis=0)
            # More efficient to use in flattened form, use indexing. Since highly sparsified, no memory issue
            # Note that fulldata is flattened
            fulldata = X.data
            left = 0
            right = 0
            for i, j in enumerate(n_nonzero_elements):
                right = right + j
                means[i] = np.mean(fulldata[left:right])
                left = right
        else:
            # non-sparse version
            n_nonzero_elements = np.count_nonzero(X, axis=0)
            means = np.mean(X, axis=0)
    else:
        # ZINB will be implemented soon
        return 0
    return means


def _tail_var_estimate(adata, mask, model='rough'):
    # Method ZINB will be implemented soon
    if model not in {'rough', 'zinb'}:
        model = 'rough'
        logg.warn('Model should be either rough or zinb (zero-inflated negative binomial)')
    X = adata.X[mask, :]
    n_cells = X.shape[0]
    n_genes = X.shape[1]
    variances = np.zeros((n_genes,))
    if model == 'rough':
        if issparse(X):
            n_nonzero_elements = X.getnnz(axis=0)
            # More efficient to use in flattened form, use indexing. Since highly sparsified, no memory issue
            # Note that fulldata is flattened
            fulldata = X.data
            left = 0
            right = 0
            for i, j in enumerate(n_nonzero_elements):
                right = right + j
                variances[i] = np.var(fulldata[left:right])
                left = right
        else:
            # non-sparse version
            n_nonzero_elements = np.count_nonzero(X, axis=0)
            variances = np.var(X, axis=0)
    else:
        # ZINB will be implemented soon
        return 0
    return variances


def _Avg_CDR(adata, mask, genes, model='rough', n_genes=None):
    # In order to get the right results, it is important to use full data,or give n_genes parameter.
    # Given an adata object and a mask (corresponding to cell selection), we want to get the CDR for
    # the cells in which certain genes are expressed. This is used to plot/visualize 'cell-volume' effects.

    if n_genes is None:
        n_genes = adata.X.shape[1]

    # Initialize list
    Summed_CDR = np.zeros((len(genes),))
    N_CDR = np.zeros((len(genes),))


    # Select nonzero-entries, genes is a list:
    adata_relevant = adata[:, genes]
    X_relevant = adata_relevant.X[mask, :]
    if issparse(X_relevant):
        indices = X_relevant.nonzero()
        # The way adata was sliced, indices should start with zero, increase one at a time
        # Get total number of expressed genes in relevant file
        N = len(indices[0])
        i = 0
        while i < N:
            Summed_CDR[indices[1][i]] += adata.X[indices[0][i], :].getnnz()
            N_CDR[[indices[1][i]]] += 1
            i = i + 1
        return (Summed_CDR / N_CDR) / n_genes
    else:
        # Non-sparse version to be implemented
        return 0

def top_ranked_group_analysis(adata, groupby, groupid, n=100, special_markers=None,
                              coloring='scores', annotate=False):
    """For one group, output a detailed chart analyzing highly ranked genes detailly.

                This is a visualization tools that helps to find significant markers and get a better understanding of
                the underlying data.

                Parameters
                ----------
                adata : :class:`~anndata.AnnData`
                    Annotated data matrix.
                groupby : `str`
                    The key of the sample grouping to consider.
                groupid: int
                    The group for which detailed analysis should be displayed.
                x : 'str'
                    x-axis labelling for plots
                y : 'str'
                    y-axis labelling for plots
                groups : `str`, `list`, optional (default: `'all'`)
                    Subset of groups, e.g. `['g1', 'g2', 'g3']`, to which comparison shall
                    be restricted. If not passed, a ranking will be generated for all
                    groups.
                n : `int`, optional (default: 100)
                    Number of datapoints in the scatterplot. If less are available, use all that are available
                special_markers: 'dict', optional (default: None)
                    If provided, this should be a dict containing a list of gene names for each group in groupby.
                    Special marked genes are highlighted in the visualization
                coloring : {'scores', 'absolute'}, optional (default: 'scores')
                    Rank either according to Scores, or to absolute test-statistic value.
                    In either case, results are scaled so as to guarantee sufficient contrast.
                annotate : bool, optional (default: False)
                    If True, annotate each datapoint in (each?) scatterplot. Helps identify specific genes. Only
                    recommended for small n.
            """
    groups = 'all'

    groups_order, groups_masks = utils.select_groups(
        adata, groups, groupby)


    imask=groupid
    mask=groups_masks[imask]

    score_list = list()
    name_list = list()
    special_markers_indices = list()
    # Note: No duplicates in each group
    for j, k in enumerate(adata.uns['rank_genes_groups_gene_scores']):
        # Make sure only first n datapoints are used
        if j >= n:
            break
        score_list.append(k[imask])
        name_list.append(adata.uns['rank_genes_groups_gene_names'][j][imask])
        # Inefficient and not generalizable but works: Check added index if in list of specially_marked_genes
        # TODO: Speed up if becomes a time issue
        if special_markers is None:
            pass
        elif adata.uns['rank_genes_groups_gene_names'][j][imask] in special_markers[imask]:
            special_markers_indices.append(len(name_list) - 1)
        else:
            pass

    ### Get all the key figures
    # make things faster by calculating only what is required for plot
    mask_rest = ~mask

    # Unlike in the above case, we need virtually all information, except for possibly CDR
    rate_group = _zero_inflation_estimate(adata[:, name_list], mask)
    rate_rest = _zero_inflation_estimate(adata[:, name_list], mask_rest)
    rate_difference= rate_group-rate_rest
    means_group = _tail_mean_estimate(adata[:, name_list], mask)
    means_rest = _tail_mean_estimate(adata[:, name_list], mask_rest)
    mean_difference=means_group-means_rest
    var_group = _tail_var_estimate(adata[:, name_list], mask)
    var_rest = _tail_var_estimate(adata[:, name_list], mask_rest)
    CDR = _Avg_CDR(adata, mask, name_list, model='rough', n_genes=None)
    adata_relevant = adata[:, name_list]
    exp, full_var_group = simple._get_mean_var(adata_relevant.X[mask])
    adata_relevant = adata[:, name_list]
    exp_rest, full_var_rest = simple._get_mean_var(adata_relevant.X[mask_rest])

    ### Prepare for coloring
    # get colored scatterplot
    # For coloring, get max score value, normalize (0,1)
    # Depending on whether normalization should be scale-invariant or only rank-invariant, do the following
    if coloring == 'scores':
        score_list = score_list / max(score_list)
        colors = cm.jet(score_list)
    elif coloring == 'absolute':
        color_list = rankdata(score_list)
        max_values = max(color_list)
        colors = cm.jet(color_list / max_values)
        # Identify true markers distinctly by using different size.
    else:
        logg.error('coloring should be either <socres> or <absolute>')
    s = 20 * np.ones(len(score_list))
    # This works for numpy access (not for normal lists though)
    s[special_markers_indices] = 100

    # Now specifically say how each subplot should look like,loop over it

    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, figsize=(12, 12))

    ax1.scatter(mean_difference, score_list, s=s, color=colors)
    ax1.set_xlabel("Expression Rate Group")
    ax1.set_ylabel("Expression Rate Rest")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax1.annotate(txt, (mean_difference[i], score_list[i]))

    ax2.scatter(rate_group, score_list, s=s, color=colors)
    ax2.set_xlabel("Expression Rate Group")
    ax2.set_ylabel("Expression Rate Rest")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax2.annotate(txt, (rate_group[i], score_list[i]))

    ax3.scatter(rate_group, rate_rest, s=s, color=colors)
    ax3.set_xlabel("Expression Rate Group")
    ax3.set_ylabel("Expression Rate Rest")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax3.annotate(txt, (rate_group[i], rate_rest[i]))

    ax4.scatter(CDR, score_list, s=s, color=colors)
    ax4.set_xlabel("Cellular Detection Rate")
    ax4.set_ylabel("Score")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax4.annotate(txt, (CDR[i], score_list[i]))

    ax5.scatter(rate_difference, score_list, s=s, color=colors)
    ax5.set_xlabel("Expression Rate Difference")
    ax5.set_ylabel("Score")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax5.annotate(txt, (rate_difference[i], score_list[i]))

    ax6.scatter(rate_group, score_list, s=s, color=colors)
    ax6.set_xlabel("Expression Rate Group")
    ax6.set_ylabel("Score")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax6.annotate(txt, (rate_group[i], score_list[i]))

    ax7.scatter(rate_group, score_list, s=s, color=colors)
    ax7.set_xlabel("Expression Rate Group")
    ax7.set_ylabel("Expression Rate Rest")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax7.annotate(txt, (rate_group[i], score_list[i]))

    ax8.scatter(rate_rest,full_var_rest, s=s, color=colors)
    ax8.set_xlabel("Expression Rate Rest")
    ax8.set_ylabel("Full Variance Rest")
    if annotate is True:
        for i, txt in enumerate(name_list):
            if i<20:
                ax8.annotate(txt, (rate_rest[i], full_var_rest[i]))
            else:
                break
    ax9.scatter(rate_group, full_var_group, s=s, color=colors)
    ax9.set_xlabel("Expression Rate Group")
    ax9.set_ylabel("Full Variance Group")
    if annotate is True:
        for i, txt in enumerate(name_list):
            ax9.annotate(txt, (rate_group[i], full_var_group[i]))

    # Here, write what we want in each subplot
    plt.tight_layout()
    plt.show()


def scatter(adata, groupby, groupid, x,y, n=100, special_markers=None,
                              coloring='scores', size=12, annotate=True):
    """For one group, output a detailed chart analyzing highly ranked genes detailly.

                This is a visualization tools that helps to find significant markers and get a better understanding of
                the underlying data.

                Parameters
                ----------
                adata : :class:`~anndata.AnnData`
                    Annotated data matrix.
                groupby : `str`
                    The key of the sample grouping to consider.
                groupid: int
                    The group for which detailed analysis should be displayed.
                x : 'str'
                    x-axis labelling for plots
                y : 'str'
                    y-axis labelling for plots
                n : `int`, optional (default: 100)
                    Number of datapoints in the scatterplot. If less are available, use all that are available
                special_markers: 'dict', optional (default: None)
                    If provided, this should be a dict containing a list of gene names for each group in groupby.
                    Special marked genes are highlighted in the visualization
                coloring : {'scores', 'absolute'}, optional (default: 'scores')
                    Rank either according to Scores, or to absolute test-statistic value.
                    In either case, results are scaled so as to guarantee sufficient contrast.
                size: int, optional (default: 12)
                    Determines scatter plot size. Large scatter-plots make it easier to identify specific genes using
                    annotate=True
                annotate : bool, optional (default: False)
                    If True, annotate each datapoint in (each?) scatterplot. Helps identify specific genes. Only
                    recommended for small n.
            """
    groups = 'all'

    groups_order, groups_masks = utils.select_groups(
        adata, groups, groupby)


    imask=groupid
    mask=groups_masks[imask]

    score_list = list()
    name_list = list()
    special_markers_indices = list()
    # Note: No duplicates in each group
    for j, k in enumerate(adata.uns['rank_genes_groups_gene_scores']):
        # Make sure only first n datapoints are used
        if j >= n:
            break
        score_list.append(k[imask])
        name_list.append(adata.uns['rank_genes_groups_gene_names'][j][imask])
        # Inefficient and not generalizable but works: Check added index if in list of specially_marked_genes
        # TODO: Speed up if becomes a time issue
        if special_markers is None:
            pass
        elif adata.uns['rank_genes_groups_gene_names'][j][imask] in special_markers[imask]:
            special_markers_indices.append(len(name_list) - 1)
        else:
            pass

    ### Get all the key figures
    # make things faster by calculating only what is required for plot
    mask_rest = ~mask

    # Get rate of expression
    rate_group = _zero_inflation_estimate(adata[:, name_list], mask)
    rate_rest = _zero_inflation_estimate(adata[:, name_list], mask_rest)
    if (x in {'full_mean_group', 'tail_mean_group', 'full_mean_difference',
              'tail_mean_difference'} or y in {'full_mean_group', 'tail_mean_group', 'full_mean_difference',
                                               'tail_mean_difference'}):
        means_group = _tail_mean_estimate(adata[:, name_list], mask)
    if (x in {'full_mean_rest', 'tail_mean_rest', 'full_mean_difference',
              'tail_mean_difference'} or y in {'full_mean_rest', 'tail_mean_rest', 'full_mean_difference',
                                               'tail_mean_difference'}):
        means_rest = _tail_mean_estimate(adata[:, name_list], mask_rest)

    if (x == 'tail_var_group' or y == 'tail_var_group'):
        # Get tail variance of expression
        var_group = _tail_var_estimate(adata[:, name_list], mask)
    if (x == 'tail_var_rest' or y == 'tail_var_rest'):
        var_rest = _tail_var_estimate(adata[:, name_list], mask_rest)
    if (x == 'CDR' or y == 'CDR'):
        # Get CDR: Need to give full adata object, since we need to count everything
        CDR = _Avg_CDR(adata, mask, name_list, model='rough', n_genes=None)
    if (x == 'full_var_group' or y == 'full_var_group'):
        # Slice first appropriately:
        adata_relevant = adata[:, name_list]
        exp, full_var_group = simple._get_mean_var(adata_relevant.X[mask])
    if (x == 'full_var_rest' or y == 'full_var_rest'):
        # Slice first appropriately:
        adata_relevant = adata[:, name_list]
        exp_rest, full_var_rest = simple._get_mean_var(adata_relevant.X[mask_rest])

    ### Prepare for coloring
    # get colored scatterplot
    # For coloring, get max score value, normalize (0,1)
    # Depending on whether normalization should be scale-invariant or only rank-invariant, do the following
    if coloring == 'scores':
        score_list = score_list / max(score_list)
        colors = cm.jet(score_list)
    elif coloring == 'absolute':
        color_list = rankdata(score_list)
        max_values = max(color_list)
        colors = cm.jet(color_list / max_values)
        # Identify true markers distinctly by using different size.
    else:
        logg.error('coloring should be either <socres> or <absolute>')
    s = 20 * np.ones(len(score_list))
    # This works for numpy access (not for normal lists though)
    s[special_markers_indices] = 100
    # In future, build method to mark top genes specially

    ### Actually do the plotting: Looping is inefficient and lengthy, but clear style
    # Potential values for x, y: 'mean' ('full' or 'tail'), 'tail_variance', 'inflation', 'CDR',
    # tail_variance_rest, Score (Just the ranking as given by test-statistic), 'full_var', 'full_var_rest'

    if x == 'expression_rate_difference':
        x_plot = rate_group - rate_rest
    elif x == 'expression_rate_group':
        x_plot = rate_group
    elif x == 'expression_rate_rest':
        x_plot = rate_rest
    elif x == 'Score':
        x_plot = score_list
    elif x == 'full_mean_difference':
        x_plot = means_group * rate_group - means_rest * rate_rest
    elif x == 'full_mean_group':
        x_plot = means_group * rate_group
    elif x == 'full_mean_rest':
        x_plot = means_rest * rate_rest
    elif x == 'tail_mean_difference':
        x_plot = means_group - means_rest
    elif x == 'tail_mean_group':
        x_plot = means_group
    elif x == 'tail_mean_rest':
        x_plot = means_rest
    elif x == 'tail_var_group':
        x_plot = var_group
    elif x == 'tail_var_rest':
        x_plot = var_rest
    elif x == 'full_var_group':
        x_plot = full_var_group
    elif x == 'full_var_rest':
        x_plot = full_var_rest
    elif x == 'CDR':
        x_plot = CDR
    else:
        logg.error('No accepted input. Check function documentation to get an overview over all inputs')

    if y == 'expression_rate_difference':
        y_plot = rate_group - rate_rest
    elif y == 'expression_rate_group':
        y_plot = rate_group
    elif y == 'expression_rate_rest':
        y_plot = rate_rest
    elif y == 'Score':
        y_plot = score_list
    elif y == 'full_mean_difference':
        y_plot = means_group * rate_group - means_rest * rate_rest
    elif y == 'full_mean_group':
        y_plot = means_group * rate_group
    elif y == 'full_mean_rest':
        y_plot = means_rest * rate_rest
    elif y == 'tail_mean_difference':
        y_plot = means_group - means_rest
    elif y == 'tail_mean_group':
        y_plot = means_group
    elif y == 'tail_mean_rest':
        y_plot = means_rest
    elif y == 'tail_var_group':
        y_plot = var_group
    elif y == 'tail_var_rest':
        y_plot = var_rest
    elif y == 'full_var_group':
        y_plot = full_var_group
    elif y == 'full_var_rest':
        y_plot =     full_var_rest
    elif y == 'CDR':
        y_plot = CDR
    else:
        logg.error('No accepted input. Check function documentation to get an overview over all inputs')

    # To make different scalings easier to compare, we set fixed limits for the case that x,y are e
    # expression rates
    if (x in {'expression_rate_difference', 'expression_rate_group', 'expression_rate_rest'} and
                y in {'expression_rate_difference', 'expression_rate_group', 'expression_rate_rest'}):
        plt.xlim(0, 1)
        plt.ylim(0, 1)
    fig, ax= plt.subplots(figsize=(size,size))
    ax.scatter(x_plot, y_plot, color=colors, s=s)
    plt.xlabel(x)
    plt.ylabel(y)
    if annotate is True:
        for i, txt in enumerate(name_list):
            plt.annotate(txt, (x_plot[i], y_plot[i]))
    plt.show()

def ROC_AUC_analysis(adata,groupby,group,n_genes=100, special_markers=None, coloring='scores', size=12, annotate=False):
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
        special_markers: 'dict', optional (default: None)
            If provided, this should be a dict containing a list of gene names for each group in groupby.
            Special marked genes are highlighted in the visualization
        coloring : {'scores', 'absolute'}, optional (default: 'scores')
            Rank either according to Scores, or to absolute test-statistic value.
            In either case, results are scaled so as to guarantee sufficient contrast.
        size: int, optional (default: 12)
            Determines scatter plot size. Large scatter-plots make it easier to identify specific genes using
            annotate=True
        annotate : bool, optional (default: False)
            If True, annotate each datapoint in (each?) scatterplot. Helps identify specific genes. Only
            recommended for small n.

    """

    fpr=adata.uns['ROCfpr' + groupby + str(group)]
    tpr=adata.uns['ROCtpr' + groupby + str(group)]
    # We dont need thresholds here

    # TODO: ALlow for comparison with rest, weighting...
    groups = 'all'

    groups_order, groups_masks = utils.select_groups(
        adata, groups, groupby)

    imask = group
    mask = groups_masks[imask]

    score_list = list()
    name_list = list()
    special_markers_indices = list()
    # Note: No duplicates in each group
    for j, k in enumerate(adata.uns['rank_genes_groups_gene_scores']):
        # Make sure only first n datapoints are used
        if j >= n_genes:
            break
        score_list.append(k[imask])
        name_list.append(adata.uns['rank_genes_groups_gene_names'][j][imask])
        # Inefficient and not generalizable but works: Check added index if in list of specially_marked_genes
        # TODO: Speed up if becomes a time issue
        if special_markers is None:
            pass
        elif adata.uns['rank_genes_groups_gene_names'][j][imask] in special_markers[imask]:
            special_markers_indices.append(len(name_list) - 1)
        else:
            pass

    ### Get all the key figures
    # make things faster by calculating only what is required for plot
    mask_rest = ~mask

    # Get rate of expression
    rate_group = _zero_inflation_estimate(adata[:, name_list], mask)
    rate_rest = _zero_inflation_estimate(adata[:, name_list], mask_rest)

    if coloring == 'scores':
        score_list = score_list / max(score_list)
        colors = cm.jet(score_list)
    elif coloring == 'absolute':
        color_list = rankdata(score_list)
        max_values = max(color_list)
        colors = cm.jet(color_list / max_values)
        # Identify true markers distinctly by using different size.
    else:
        logg.error('coloring should be either <socres> or <absolute>')
    s = 20 * np.ones(len(score_list))
    # This works for numpy access (not for normal lists though)
    s[special_markers_indices] = 100

    fig, ax = plt.subplots(figsize=(size, size))
    ax.scatter(rate_rest, rate_group, color=colors, s=s)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    if annotate is True:
        for i, txt in enumerate(name_list):
            plt.annotate(txt, (rate_rest[i], rate_group[i]))
            # TODO: Add AUC

    # Now annotate the lines for all genes:
    # TODO: Until now, everything requires same number of n (i.e. all in name list). Shouldn't be the case. Resolve.
    for i,j in enumerate(name_list):
        plt.plot(fpr[name_list[i]], tpr[name_list[i]], color=colors[i])
    plt.show()

def comparison_table(adata, name_keys, group=None, color_thresholds=None, n_genes=70):
    ## Very early: Take list of keys for which to search adata.uns.
    ## Then build a table of all included genes, see which w
    ## Trick: Map only the top 20/30 or so, color green/yellow/orange/red for top 40/ 60 / below
    ##TODO: Add functionality for group, color thresholds.

    ## Assume all annotations have the same length
    name_list={}
    for i,j in enumerate(name_keys):
        name_list[i]=adata.uns[j]

    # We dont need rank list as we assume that name list is ordered.

    length=len(name_list[0])
    width=len(name_list)
    # we create one large table (still no memory issue (realistic max. 10000*10--> 1 Mbyte approx). Truncuate later
    rank_table=length* np.ones((length*width,width))
    # Create full name list
    full_name_list=list()
    n=-1

    for key in name_list:
        for k, l in enumerate(name_list[key]):
            # Only plot for each group the top 20 genes to avoid that the table becomes to large. Max_n should be
            # a parameter in the future
            # Problem: We should add everything, but only plot a certain truncu
            if l not in full_name_list:
                full_name_list.append(l)
                n=n+1
                m=n
            else:
                m=full_name_list.index(l)
            rank_table[m,key]=k

    # Create table with all entries
    if max_n<n:
        n=max_n
    trunc_table=rank_table[0:n+1,:]
    # Do the colorings:
    colors=trunc_table.copy()
    for i in range(n+1):
        # Here, for now we use the convention that the respective levels are true if the minimum rank is larger
        # than the specified number
        top20=True
        top50=True
        top100=True
        for j in range(width):
            if colors[i,j]>=100:
                top20=top50=top100=False
            elif colors[i,j]>=51:
                top20=top50=False
            elif colors[i,j]>=21:
                top20=False
        # Now depending on the boolean values, define colors.
        if top100 is False:
            colors[i,:]=0
        elif top50 is False:
            colors[i,:]=0.5
        elif top20 is False:
            colors[i,:]=0.8
        else:
            colors[i,:]=1
    fig,ax=plt.subplots(1,1)
    ax.table(cellText=trunc_table, rowLabels=full_name_list[0:n+1], colLabels=name_keys,
                        cellColours=cm.brg(colors))
    plt.tight_layout()
    plt.show()




def comparison_v2(adata, name_keys, group=None, color_thresholds=None, n_genes=70):
    name_list_cut = {}
    for i, j in enumerate(name_keys):
        name_list_cut[i] = adata.uns[j][0:n_genes]
    name_list = {}
    for i, j in enumerate(name_keys):
        name_list[i] = adata.uns[j]

    length = n_genes
    width = len(name_list)

    rank_table = pd.DataFrame(name_list_cut)
    row_names=np.arange(n_genes)+1
    colors=np.ndarray((length,width))
    for key in name_list:
        for i in range(n_genes):
            top100=False
            top50=False
            top20=False
            for key2 in name_list:
                if key is key2:
                    pass
                else:
                    if name_list[key][i] in name_list[key2]:
                        index=name_list[key2].index(name_list[key][i])
                        if index <100:
                            top100=True
                            if index <50:
                                top50=True
                                if index >=20:
                                    top20=True
                        else:
                            pass
                    else:
                        top100=False
                        top50=False
                        top20=False
            if top100 is True:
                colors[i, key] = 0.55
                if top50 is True:
                    colors[i, key] = 0.75
                    if top20 is True:
                        colors[i, key] = 0.9
            else:
                colors[i, :] = 0.35

    plt.figure(figsize=(4,4 ), dpi=120)
    ax = plt.subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.table(cellText=rank_table.as_matrix(), rowLabels=row_names, colLabels=name_keys,
             cellColours=cm.afmhot(colors), loc="center", fontsize=22)
    plt.show()
