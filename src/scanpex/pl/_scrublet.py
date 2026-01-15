import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc

from .. import pp


def scrublet(
    adata: ad.AnnData,
    ax: plt.Axes,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    **kwargs
):
    """
    Visualize predicted doublets on quality control metrics.

    This function checks if the Scrublet prediction exists in `adata.obs`.
    If not, it executes the prediction using `pp.scrublet`. It then generates
    a scatter plot where predicted doublets are highlighted. The data is sorted
    prior to plotting to ensure that doublets (True) are plotted on top of
    singlets (False) for better visibility.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix.
    ax : matplotlib.axes.Axes
        The axis on which to draw the scatter plot.
    x : str, optional
        The column name in `adata.obs` for the x-axis.
        By default "total_counts".
    y : str, optional
        The column name in `adata.obs` for the y-axis.
        By default "n_genes_by_counts".
    **kwargs
        Additional keyword arguments passed to `scanpy.pl.scatter`.
        If `palette` is not provided, it defaults to coloring singlets
        "lightgrey" and doublets "tab:red".

    Returns
    -------
    None
        The plot is drawn directly onto the provided `ax` object.
    """
    if "predicted_doublet" not in adata.obs.columns:
        pp.scrublet(adata=adata, remove=False, update=False)

    order = adata.obs["predicted_doublet"].argsort()
    palette = kwargs.pop("palette", {False: "lightgrey", True: "tab:red"})

    sc.pl.scatter(
        adata[order],
        x=x,
        y=y,
        color="predicted_doublet",
        show=False,
        ax=ax,
        palette=palette,
        **kwargs
    )
