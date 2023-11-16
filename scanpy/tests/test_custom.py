from test_highly_variable_genes import (
    test_highly_variable_genes_subset_inplace_consistency,
)


for flavor in [
    "seurat",
]:  # "cell_ranger", "seurat_v3"]:
    for subset in [True, False]:
        for inplace in [True, False]:
            test_highly_variable_genes_subset_inplace_consistency(
                flavor, subset, inplace
            )


a = 0
