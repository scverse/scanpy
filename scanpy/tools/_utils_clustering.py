def rename_groups(adata, key_added, restrict_key, restrict_categories,
        restrict_indices, groups):
    key_added = restrict_key + '_R' if key_added is None else key_added
    all_groups = adata.obs[restrict_key].astype('U')
    prefix = '-'.join(restrict_categories) + ','
    new_groups = [prefix + g for g in groups.astype('U')]
    all_groups.iloc[restrict_indices] = new_groups
    return all_groups


def restrict_adjacency(adata, restrict_key, restrict_categories, adjacency):
    if not isinstance(restrict_categories[0], str):
        raise ValueError('You need to use strings to label categories, '
                         'e.g. \'1\' instead of 1.')
    for c in restrict_categories:
        if c not in adata.obs[restrict_key].cat.categories:
            raise ValueError(
                '\'{}\' is not a valid category for \'{}\''
                .format(c, restrict_key))
    restrict_indices = adata.obs[restrict_key].isin(restrict_categories).values
    adjacency = adjacency[restrict_indices, :]
    adjacency = adjacency[:, restrict_indices]
    return adjacency, restrict_indices

