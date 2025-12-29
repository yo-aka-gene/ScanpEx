from typing import Optional

import anndata as ad
import scanpy as sc


def scrublet(
    adata: ad.AnnData, remove: bool = False, update: bool = False, **kwargs
) -> Optional[ad.AnnData]:
    if ("predicted_doublet" not in adata.obs.columns) or update:
        sc.pp.scrublet(adata, **kwargs)

    if remove:
        return adata[~adata.obs["predicted_doublet"], :].copy()

    return None
