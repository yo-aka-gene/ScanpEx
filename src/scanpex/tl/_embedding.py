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
        adata_return.uns["neighbors"] = copy.deepcopy(adata_source.uns["neighbors"])

    for key, value in adata_source.obsp.items():
        adata_return.obsp[key] = value.copy()

    return None if overwrite else adata_return


def transfer_clustering_info(
    adata_source: ad.AnnData,
    adata_target: ad.AnnData,
    obs_keys: str | list[str] = "leiden",
    overwrite: bool = False,
) -> ad.AnnData | None:
    """
    Transfer clustering annotations and unstructured metadata between AnnData objects.

    This function copies selected observation-level annotations from a source
    `AnnData` object to a target `AnnData` object together with all entries in
    `uns`.

    It is intended for workflows where clustering or downstream annotation is
    performed on a derived `AnnData` object, while the resulting annotations
    need to be transferred back to the original object.

    Because observation-level annotations are indexed by observations, the
    source and target objects must have identical `obs_names` in the same
    order.

    Parameters
    ----------
    adata_source : anndata.AnnData
        Source AnnData object containing the clustering annotations and
        unstructured metadata to transfer.
    adata_target : anndata.AnnData
        Target AnnData object receiving the transferred information. Its
        observations must match `adata_source.obs_names` exactly and in the
        same order.
    obs_keys : str or list of str, default "leiden"
        Observation metadata column or columns in `adata_source.obs` to
        transfer to `adata_target.obs`.
    overwrite : bool, default False
        Whether to modify `adata_target` in place. If False, a copy of
        `adata_target` is created and returned. If True, `adata_target` is
        modified directly and the function returns None.

    Returns
    -------
    anndata.AnnData or None
        A copy of `adata_target` containing the transferred annotations and
        unstructured metadata when `overwrite=False`. Returns None when
        `overwrite=True`.

    Raises
    ------
    ValueError
        If `adata_source` and `adata_target` do not have identical
        `obs_names` in the same order.
    KeyError
        If any requested key is not present in `adata_source.obs`.

    Notes
    -----
    Observation metadata columns are copied individually, while `uns` is
    deep-copied in its entirety to avoid sharing mutable objects between the
    source and target AnnData objects.
    """
    if not adata_source.obs_names.equals(adata_target.obs_names):
        raise ValueError(
            "`adata_source` and `adata_target` must have identical "
            "`obs_names` in the same order."
        )

    if isinstance(obs_keys, str):
        obs_keys = [obs_keys]

    missing_keys = [key for key in obs_keys if key not in adata_source.obs.columns]
    if missing_keys:
        raise KeyError(
            f"The following keys are not present in `adata_source.obs`: "
            f"{missing_keys}"
        )

    adata_return = adata_target if overwrite else adata_target.copy()

    for key in obs_keys:
        adata_return.obs[key] = adata_source.obs[key].copy()

    adata_return.uns = copy.deepcopy(adata_source.uns)

    return None if overwrite else adata_return
