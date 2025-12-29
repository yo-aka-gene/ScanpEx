from typing import Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import zscore


def sigmoid(x: Union[np.ndarray, pd.Series]) -> np.ndarray:
    """
    Apply sigmoid transformation subsequently to the Z-score transformation.

    Note that Z-score calculation requires a distribution (array-like),
    not a single scalar value.

    Parameters
    ----------
    x : Union[np.ndarray, pd.Series]
        Array-like data to be converted (e.g., a column of adata.obs).

    Returns
    -------
    np.ndarray
        Sigmoid-transformed Z-score of x (values between 0 and 1).
    """
    return 1 / (1 + np.exp(-zscore(x, nan_policy="omit")))


def prob_genes(adata: ad.AnnData, gene_list: list, **kwargs) -> Union[None, ad.AnnData]:
    """
    Calculate gene signature scores and transform them to [0, 1] probability
    using the sigmoid-Z-score transformation.
    This is a wrapper of `scanpy.tl.score_genes`.
    It adds a `{score_name}_prob` column to `adata.obs`.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix.
    gene_list : list
        The list of genes to be scored.
    **kwargs
        Additional arguments passed to `scanpy.tl.score_genes`.
        (e.g., `score_name`, `ctrl_size`, `random_state`)

    Returns
    -------
    Union[None, ad.AnnData]
        Returns `None` if `copy=False` (default), otherwise returns a copy of `adata`.
        The result is stored in `adata.obs['{score_name}_prob']`.
    """
    score_name = kwargs.get("score_name", "score")
    result = sc.tl.score_genes(adata=adata, gene_list=gene_list, **kwargs)

    data = result if result is not None else adata
    data.obs[f"{score_name}_prob"] = sigmoid(data.obs[score_name])

    return result


def score_genes_cell_cycle(
    adata: ad.AnnData, s_genes: list, g2m_genes: list, **kwargs
) -> None:
    """
    Assign cell cycle phases based on "S_score" and "G2M_score".

    This is a wrapper of `scanpex.tl.prob_genes` and `scanpy.tl.score_genes_cell_cycle`.
    It calculates raw scores, probability scores (0-1), and assigns phases.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix.
    s_genes : list
        The list of S-phase related genes to be scored.
    g2m_genes : list
        The list of G2/M-phase related genes to be scored.
    **kwargs
        Additional arguments passed to `scanpy.tl.score_genes` and
        `scanpy.tl.score_genes_cell_cycle`.

    Returns
    -------
    None
        The following columns are added to `adata.obs`:
        - "S_score", "S_score_prob"
        - "G2M_score", "G2M_score_prob"
        - "phase" (with "G2M" label prettified to "G2/M")
    """
    prob_genes(
        adata=adata, gene_list=s_genes, score_name="S_score", copy=False, **kwargs
    )
    prob_genes(
        adata=adata, gene_list=g2m_genes, score_name="G2M_score", copy=False, **kwargs
    )
    sc.tl.score_genes_cell_cycle(
        adata=adata, s_genes=s_genes, g2m_genes=g2m_genes, copy=False, **kwargs
    )
    adata.obs["phase"] = adata.obs["phase"].map(lambda p: p if p != "G2M" else "G2/M")


def curate_phase(adata: ad.AnnData, thresh: float = 0) -> None:
    """
    Manually curate the threshold for cell cycle phase assignment.

    This function re-assigns phases (G1, S, G2/M) based on existing "S_score"
    and "G2M_score" using a user-defined threshold.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with pre-assigned cell cycle scores
        (e.g., via `scanpex.tl.score_genes_cell_cycle`).
    thresh : float, optional (default: 0)
        The threshold value for determination.
        Scores below this value are considered G1 phase.

    Returns
    -------
    None
        `adata.obs["phase"]` is overwritten with the new assignments.
    """
    # Create flags using decimal place value logic
    # 10s place: S score exceeds threshold
    s_phase = adata.obs["S_score"].map(lambda s: 10 if s >= thresh else 0)

    # 100s place: G2M score exceeds threshold
    g2m_phase = adata.obs["G2M_score"].map(lambda s: 100 if s >= thresh else 0)

    # 1s place: G2M score is strictly larger than S score
    g2m_is_larger = (adata.obs["G2M_score"] - adata.obs["S_score"]).map(
        lambda s: 1 if s > 0 else 0
    )

    # Calculate unique key for each state
    # 0   (000) -> G1 (Both < thresh)
    # 1   (001) -> G1 (Both < thresh, G2M > S) ... unlikely but logically G1
    # 10  (010) -> S  (S >= thresh > G2M)
    # 101 (101) -> G2/M (G2M >= thresh > S)
    # 110 (110) -> S    (Both >= thresh, S >= G2M)
    # 111 (111) -> G2/M (Both >= thresh, G2M > S)
    adata.obs["phase"] = (
        (s_phase + g2m_phase + g2m_is_larger)
        .map(
            lambda s: {0: "G1", 1: "G1", 10: "S", 101: "G2/M", 110: "S", 111: "G2/M"}[s]
        )
        .astype("category")
    )
