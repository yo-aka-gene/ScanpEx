from typing import Union

import anndata as ad
import numpy as np
import scanpy as sc
from scipy.stats import zscore


def prob_genes(
    adata: ad.AnnData,
    gene_list: list,
    **kwargs
) -> Union[None, ad.AnnData]:
    score_name = kwargs.get("score_name", "score")
    result = sc.tl.score_genes(
        adata=adata, 
        gene_list=gene_list,
        **kwargs
    )

    sigmoid = lambda s: 1 / (1 + np.exp(-zscore(s, nan_policy="omit")))

    data = result if result is not None else adata
    data.obs[f"{score_name}_prob"] = sigmoid(data.obs[score_name])

    return result


def score_genes_cell_cycle(
    adata: ad.AnnData,
    s_genes: list,
    g2m_genes: list,
    **kwargs
) -> None:
    prob_genes(
        adata=adata, 
        gene_list=s_genes, 
        score_name='S_score', 
        copy=False,
        **kwargs
    )
    prob_genes(
        adata=adata, 
        gene_list=g2m_genes, 
        score_name='G2M_score', 
        copy=False,
        **kwargs
    )
    sc.tl.score_genes_cell_cycle(
        adata=adata,
        s_genes=s_genes, 
        g2m_genes=g2m_genes,
        copy=False,
        **kwargs
    )
    adata.obs["phase"] = adata.obs["phase"].map(lambda p: p if p != "G2M" else "G2/M")


def currate_phase(
        adata: ad.AnnData, 
        thresh: float = 0
    ) -> None:
    s_phase = adata.obs["S_score"].map(lambda s: 10 if s >= thresh else 0)
    g2m_phase = adata.obs["G2M_score"].map(lambda s: 100 if s >= thresh else 0)
    g2m_is_larger = (adata.obs["G2M_score"] - adata.obs["S_score"]).map(
        lambda s: 1 if s > 0 else 0
    )
    adata.obs["phase"] = (s_phase + g2m_phase + g2m_is_larger).map(
        lambda s: {0: "G1", 1: "G1", 10: "S", 101: "G2/M", 110: "S", 111: "G2/M"}[s]
    )
