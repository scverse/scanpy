Remove `scanpy.pp.recipe_seurat`, `scanpy.pp.recipe_zheng17`, and `scanpy.pp.recipe_weinreb17` and their private implementation.
Use explicit preprocessing steps instead; see {doc}`/how-to/legacy-recipes` for migration examples and the numerical differences to consider for Weinreb workflows.
{func}`scanpy.experimental.pp.recipe_pearson_residuals` remains available.
