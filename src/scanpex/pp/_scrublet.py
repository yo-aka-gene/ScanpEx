from typing import Optional

import anndata as ad
import scanpy as sc


def scrublet(
    adata: ad.AnnData, remove: bool = False, update: bool = False, **kwargs
) -> Optional[ad.AnnData]:
    """
    Run Scrublet to predict and optionally remove doublets.

    This function wraps `scanpy.pp.scrublet`. It checks if doublet prediction
    has already been performed to avoid redundant computation unless `update`
    is set to True.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix of shape (n_obs, n_vars).
    remove : bool, optional
        If True, returns a new AnnData object with predicted doublets removed.
        If False, modifies `adata` in-place by adding 'predicted_doublet'
        column to `.obs`. By default False.
    update : bool, optional
        If True, forces re-execution of Scrublet even if 'predicted_doublet'
        already exists in `adata.obs`. By default False.
    **kwargs
        Additional keyword arguments passed to `scanpy.pp.scrublet`.

    Returns
    -------
    ad.AnnData or None
        If `remove` is True, returns a subsetted AnnData object with doublets
        filtered out. Otherwise, returns None (updates `adata` in-place).
    """
    if ("predicted_doublet" not in adata.obs.columns) or update:
        sc.pp.scrublet(adata, **kwargs)

    if remove:
        return adata[~adata.obs["predicted_doublet"], :].copy()

    return None
