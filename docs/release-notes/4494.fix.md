Make {func}`scanpy.pl.scatter` look up `basis` in `adata.obsm` exactly before falling back to an `X_` prefix, consistent with {func}`scanpy.pl.embedding`.
Support custom and explicitly prefixed keys while preserving diffusion component numbering for `diffmap` and `X_diffmap`.
