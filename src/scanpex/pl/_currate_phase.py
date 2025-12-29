import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from scanpex.tl import currate_phase as cp


def currate_phase(
    adata: ad.AnnData,
    thresh: float,
    ax: plt.Axes,
    s: int = 1,
    lw: float = .5,
    color: str = ".2"
) -> None:
    cp(adata=adata, thresh=thresh)
    sns.scatterplot(data=adata.obs, x="S_score", y="G2M_score", hue="phase", ax=ax, s=s)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    ax.vlines(thresh, ylim[0], thresh, color=color, lw=lw)
    ax.hlines(thresh, xlim[0], thresh, color=color, lw=lw)
    ax.plot(
        [thresh, max(xlim[1], ylim[1])], [thresh, max(xlim[1], ylim[1])],
        color=color, lw=lw
    )
    ax.set_xlim(*xlim), ax.set_ylim(*ylim)
    ax.legend(loc="center left", bbox_to_anchor=(1, .5), title="Phase", frameon=False)
    ax.set_aspect('equal');
