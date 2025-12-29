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
