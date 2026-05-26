import anndata as ad
import numpy as np


def get_degs(
    data: ad.AnnData,
    cluster_name: str,
    pvals_adj: float = 0.05,
    max_logfc: float = np.inf,
    min_logfc: float = -np.inf,
) -> np.ndarray:
    loc = (
        (data.uns["rank_genes_groups"]["pvals_adj"][cluster_name] < pvals_adj)
        & (data.uns["rank_genes_groups"]["logfoldchanges"][cluster_name] < max_logfc)
        & (data.uns["rank_genes_groups"]["logfoldchanges"][cluster_name] > min_logfc)
    )
    return data.uns["rank_genes_groups"]["names"][cluster_name][loc]


def top_degs(
    data: ad.AnnData,
    top: int = 5,
    pvals_adj: float = 0.05,
    max_logfc: float = np.inf,
    min_logfc: float = -np.inf,
    unique: bool = True,
) -> np.ndarray:
    arr = np.concatenate(
        [
            get_degs(
                data=data,
                cluster_name=name,
                pvals_adj=pvals_adj,
                max_logfc=max_logfc,
                min_logfc=min_logfc,
            )[:top]
            for name in data.uns["rank_genes_groups"]["names"].dtype.names
        ],
        axis=0,
    )
    _, idx = np.unique(arr, return_index=True)
    return arr[np.sort(idx)] if unique else arr
