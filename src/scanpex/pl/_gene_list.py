from typing import Dict, Optional, Tuple, Union

import seaborn as sns

from .. import ft


def gene_list(
    gl: ft.GeneList,
    group_key: Optional[str] = "SEACells",
    method: str = "pearson",
    abs_corr: bool = False,
    use_gene_names: bool = True, 
    cbar_pos: Tuple[float, float, float, float] = (1, .4, .025, .2), 
    cbar_kws: Optional[Dict[str, str]] = None,
    dendrogram_ratio: Union[float, Tuple[float, float]] = 0.1, 
    cmap: str = "RdYlBu_r",
    rasterized: bool = True,
    **kwargs
):
    data = gl._get_aggregated_df(group_key=group_key, score_name="")
    data = data.drop(columns=["score"]) if "score" in data.columns else data
    data.columns.name = ""
    ticklabels = gl.genes if use_gene_names else gl.ids
    if cbar_kws is None:
        cbar_kws = {"label": r"$|\rho|$" if abs_corr else r"$\rho$"}

    sns.clustermap(
        data.corr(method=method).abs() if abs_corr else data.corr(method=method),
        xticklabels=ticklabels, yticklabels=ticklabels,
        cbar_pos=cbar_pos, cbar_kws=cbar_kws, dendrogram_ratio=dendrogram_ratio,
        cmap=cmap, rasterized=rasterized, **kwargs
    )
