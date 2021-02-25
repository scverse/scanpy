import numpy as np
from scanpy.tools import rank_genes_groups, filter_rank_genes_groups
from scanpy.datasets import pbmc68k_reduced


names_no_reference = np.array(
    [
        ['CD3D', 'ITM2A', 'CD3D', 'CCL5', 'CD7', 'nan', 'CD79A', 'nan', 'NKG7', 'LYZ'],
        ['CD3E', 'CD3D', 'nan', 'NKG7', 'CD3D', 'AIF1', 'CD79B', 'nan', 'GNLY', 'CST3'],
        ['IL32', 'RPL39', 'nan', 'CST7', 'nan', 'nan', 'nan', 'SNHG7', 'CD7', 'nan'],
        ['nan', 'SRSF7', 'IL32', 'GZMA', 'nan', 'LST1', 'IGJ', 'nan', 'CTSW', 'nan'],
        [
            'nan',
            'nan',
            'CD2',
            'CTSW',
            'CD8B',
            'TYROBP',
            'ISG20',
            'SNHG8',
            'GZMB',
            'nan',
        ],
    ]
)

names_reference = np.array(
    [
        ['CD3D', 'ITM2A', 'CD3D', 'nan', 'CD3D', 'nan', 'CD79A', 'nan', 'CD7'],
        ['nan', 'nan', 'nan', 'CD3D', 'nan', 'AIF1', 'nan', 'nan', 'NKG7'],
        ['nan', 'nan', 'nan', 'NKG7', 'nan', 'FCGR3A', 'ISG20', 'SNHG7', 'CTSW'],
        ['nan', 'CD3D', 'nan', 'CCL5', 'CD7', 'nan', 'CD79B', 'nan', 'GNLY'],
        ['CD3E', 'IL32', 'nan', 'IL32', 'CD27', 'FCER1G', 'nan', 'nan', 'nan'],
    ]
)


def test_filter_rank_genes_groups():
    adata = pbmc68k_reduced()

    # fix filter defaults
    args = {
        'adata': adata,
        'key_added': 'rank_genes_groups_filtered',
        'min_in_group_fraction': 0.25,
        'min_fold_change': 1,
        'max_out_group_fraction': 0.5,
    }

    rank_genes_groups(
        adata, 'bulk_labels', reference='Dendritic', method='wilcoxon', n_genes=5
    )
    filter_rank_genes_groups(**args)

    assert np.array_equal(
        names_reference,
        np.array(adata.uns['rank_genes_groups_filtered']['names'].tolist()),
    )

    rank_genes_groups(adata, 'bulk_labels', method='wilcoxon', n_genes=5)
    filter_rank_genes_groups(**args)

    assert np.array_equal(
        names_no_reference,
        np.array(adata.uns['rank_genes_groups_filtered']['names'].tolist()),
    )

    rank_genes_groups(adata, 'bulk_labels', method='wilcoxon', pts=True, n_genes=5)
    filter_rank_genes_groups(**args)

    assert np.array_equal(
        names_no_reference,
        np.array(adata.uns['rank_genes_groups_filtered']['names'].tolist()),
    )
