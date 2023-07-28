"""\
Cell-state density computation using Mellon
"""

from typing import Optional
from anndata import AnnData
import numpy as np

from ... import logging as logg


def mellon(
    adata: AnnData,
    repr_key: str = 'X_palantir_diff_comp',
    density_key: str = 'mellon_log_density',
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """\
    Compute cell-state density with Mellon.

    This function uses Mellon to compute the density of cell states,
    which is stored in the obs attribute of the AnnData object. The function returns
    the computed density. If `repr_key` is not found in the AnnData object,
    an error is raised suggesting the user to run the function `scanpy.external.tl.palantir(adata)`.

    Additionally, the density prediction model is serialized and stored in the `.uns` attribute
    of the AnnData object under the key `'{density_key}_predictor'`. This can be deserialized
    and used for prediction by using the `palantir.utils.run_density_evaluation` method.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    repr_key
        Key to retrieve cell-state representation from the AnnData object.
        Default is 'X_palantir_diff_comp'.
    density_key
        Key under which the computed density values are stored in the obs
        of the AnnData object. Default is 'mellon_log_density'.
    copy
        If `True`, return a copy of `adata` instead of writing to it.
    **kwargs
        Additional keyword arguments to be passed to mellon.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields:

    adata.obs[density_key]
        Array of log density values computed for each cell.
    adata.obs[density_key + "_clipped"]
        Array of log density values clipped between the 1st and 100th percentiles.
    adata.uns[density_key + "_predictor"]
        Serialized density prediction model.
    """
    _check_import()
    from mellon import DensityEstimator

    adata = adata.copy() if copy else adata

    logg.info('Mellon Density Maps in progress ...')

    # Check if 'X_palantir_diff_comp' key is in the adata object
    if repr_key not in adata.obsm.keys():
        raise ValueError(
            f"{repr_key} not found in adata.obsm. "
            "Please run `scanpy.external.tl.palantir(adata)` first."
        )

    X = adata.obsm[repr_key]

    # Set the default arguments for mellon.DensityEstimator
    mellon_args = dict()
    mellon_args.update(kwargs)

    dest = DensityEstimator(**mellon_args)
    log_density = np.asarray(dest.fit_predict(X))

    # Save the density values and predictor to adata
    adata.obs[density_key] = log_density
    adata.obs[density_key + "_clipped"] = np.clip(
        log_density, *np.quantile(log_density, [0.01, 1])
    )
    adata.uns[density_key + "_predictor"] = dest.predict.to_dict()

    return adata if copy else None


def _check_import():
    try:
        import mellon
    except ImportError:
        raise ImportError('\nplease install mellon:\n\tpip install mellon')
