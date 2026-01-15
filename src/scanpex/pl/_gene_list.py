from typing import Dict, Optional, Tuple, Union

import seaborn as sns

from .. import ft


def gene_list(
    gl: ft.GeneList,
    group_key: Optional[str] = "leiden",
    method: str = "pearson",
    abs_corr: bool = False,
    use_gene_names: bool = True,
    cbar_pos: Tuple[float, float, float, float] = (1, 0.4, 0.025, 0.2),
    cbar_kws: Optional[Dict[str, str]] = None,
    dendrogram_ratio: Union[float, Tuple[float, float]] = 0.1,
    cmap: str = "RdYlBu_r",
    rasterized: bool = True,
    **kwargs
):
    """
    Compute and plot the correlation matrix of genes in the GeneList.

    This function aggregates expression data based on the provided `group_key`,
    calculates pairwise correlations between genes using the specified `method`,
    and visualizes the result as a hierarchically clustered heatmap using
    `seaborn.clustermap`.

    Parameters
    ----------
    gl : ft.GeneList
        The GeneList object containing the genes of interest and methods to
        retrieve aggregated data.
    group_key : str, optional
        The key in the observation metadata used to aggregate the data
        (e.g., cell clusters). Passed to `gl._get_aggregated_df`.
        By default "leiden".
    method : str, optional
        The method to compute correlation (e.g., 'pearson', 'spearman', 'kendall').
        Passed to `pandas.DataFrame.corr`. By default "pearson".
    abs_corr : bool, optional
        If True, plots the absolute value of the correlation coefficients.
        Useful for focusing on the strength of the relationship regardless of direction.
        By default False.
    use_gene_names : bool, optional
        If True, uses gene symbols (`gl.genes`) as axis labels.
        If False, uses gene IDs (`gl.ids`). By default True.
    cbar_pos : tuple of float, optional
        The position of the colorbar (left, bottom, width, height).
        By default (1, 0.4, 0.025, 0.2).
    cbar_kws : dict, optional
        Keyword arguments for the colorbar. If None, automatically sets the
        label to "rho" or "|rho|".
    dendrogram_ratio : float or tuple of float, optional
        Proportion of the figure size devoted to the dendrograms.
        By default 0.1.
    cmap : str, optional
        The mapping from data values to color space. By default "RdYlBu_r".
    rasterized : bool, optional
        If True, rasterizes the heatmap mesh to reduce file size when saving
        as vector graphics (PDF/SVG). By default True.
    kwargs : dict
        Additional keyword arguments passed to `seaborn.clustermap`.

    Returns
    -------
    seaborn.matrix.ClusterGrid
        The ClusterGrid object returned by `seaborn.clustermap`.
    """
    data = gl._get_aggregated_df(group_key=group_key, score_name="")
    data = data.drop(columns=["score"]) if "score" in data.columns else data
    data.columns.name = ""
    ticklabels = gl.genes if use_gene_names else gl.ids
    if cbar_kws is None:
        cbar_kws = {"label": r"$|\rho|$" if abs_corr else r"$\rho$"}

    return sns.clustermap(
        data.corr(method=method).abs() if abs_corr else data.corr(method=method),
        xticklabels=ticklabels,
        yticklabels=ticklabels,
        cbar_pos=cbar_pos,
        cbar_kws=cbar_kws,
        dendrogram_ratio=dendrogram_ratio,
        cmap=cmap,
        rasterized=rasterized,
        **kwargs
    )
