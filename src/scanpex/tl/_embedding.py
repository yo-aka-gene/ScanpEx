import copy

import anndata as ad


def transfer_embedding_info(
    adata_source: ad.AnnData,
    adata_target: ad.AnnData,
    overwrite: bool = False,
) -> ad.AnnData | None:
    """Transfer embedding and neighborhood information between AnnData objects."""
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
