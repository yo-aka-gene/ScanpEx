import copy

import anndata as ad


def transfer_embedding_info(
    adata_source: ad.AnnData,
    adata_target: ad.AnnData,
    overwrite: bool = False,
) -> ad.AnnData | None:
    """
    Transfer embedding and neighborhood information between AnnData objects.

    This function copies commonly used low-dimensional representations and
    neighborhood graph information from a source `AnnData` object to a target
    `AnnData` object. It is intended for cases where multiple `AnnData` objects
    contain the same observations in the same order but differ in their
    expression matrix, annotations, or other metadata.

    The following information is transferred:

    - `obsm["X_pca"]`
    - `obsm["X_umap"]`
    - `uns["neighbors"]`, if present
    - all entries in `obsp`

    Because embeddings and neighborhood graphs are indexed by observations,
    the source and target objects must have identical `obs_names` in the same
    order.

    Parameters
    ----------
    adata_source : anndata.AnnData
        Source AnnData object containing the embedding and neighborhood
        information to transfer.
    adata_target : anndata.AnnData
        Target AnnData object receiving the transferred information. Its
        observations must match `adata_source.obs_names` exactly and in the
        same order.
    overwrite : bool, default False
        Whether to modify `adata_target` in place. If False, a copy of
        `adata_target` is created and returned. If True, `adata_target` is
        modified directly and the function returns None.

    Returns
    -------
    anndata.AnnData or None
        A copy of `adata_target` containing the transferred embedding and
        neighborhood information when `overwrite=False`. Returns None when
        `overwrite=True`.

    Raises
    ------
    ValueError
        If `adata_source` and `adata_target` do not have identical
        `obs_names` in the same order.
    KeyError
        If `adata_source` does not contain `obsm["X_pca"]` or
        `obsm["X_umap"]`.

    Notes
    -----
    The PCA and UMAP embeddings and all `obsp` entries are copied to avoid
    sharing mutable array or sparse-matrix objects between the source and
    target AnnData objects. The `uns["neighbors"]` entry is deep-copied because
    it may contain nested mutable objects.

    This function does not transfer the expression matrix, observation
    metadata, variable metadata, or arbitrary entries from `uns` or `obsm`
    other than those listed above.
    """
    if not adata_source.obs_names.equals(adata_target.obs_names):
        raise ValueError(
            "`adata_source` and `adata_target` must have identical "
            "`obs_names` in the same order."
        )

    adata_return = adata_target if overwrite else adata_target.copy()

    adata_return.obsm["X_pca"] = adata_source.obsm["X_pca"].copy()
    adata_return.obsm["X_umap"] = adata_source.obsm["X_umap"].copy()

    if "neighbors" in adata_source.uns:
        adata_return.uns["neighbors"] = copy.deepcopy(
            adata_source.uns["neighbors"]
        )

    for key, value in adata_source.obsp.items():
        adata_return.obsp[key] = value.copy()

    return None if overwrite else adata_return
