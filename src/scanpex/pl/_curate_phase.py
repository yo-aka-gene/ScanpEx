import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from scanpex.tl import curate_phase as cp


def curate_phase(
    adata: ad.AnnData,
    thresh: float,
    ax: plt.Axes,
    s: int = 1,
    lw: float = 0.5,
    color: str = ".2",
) -> None:
    """
    Assign cell cycle phases and visualize the classification boundaries.

    This function first executes the phase curation logic (via `scanpex.tl.curate_phase`)
    using the specified threshold. It then generates a scatter plot of 'S_score'
    versus 'G2M_score' and overlays the linear boundaries that demarcate the
    G1, S, and G2M phases.

    The decision boundaries are drawn as follows:
    - Vertical line: Separates G1/G2M from S (at x = `thresh`).
    - Horizontal line: Separates G1/S from G2M (at y = `thresh`).
    - Diagonal line: Separates S from G2M when both scores are high.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. Must contain 'S_score' and 'G2M_score'
        in `.obs`. The object is modified in-place to add the 'phase' column.
    thresh : float
        The score threshold used to distinguish proliferating cells from G1.
        This value determines the position of the boundary lines.
    ax : matplotlib.axes.Axes
        The axis on which to draw the plot.
    s : int, optional
        The marker size for the scatter plot points. By default 1.
    lw : float, optional
        The line width for the decision boundary lines. By default 0.5.
    color : str, optional
        The color of the decision boundary lines. By default ".2" (dark grey).

    Returns
    -------
    None
        The plot is drawn directly onto the provided `ax` object.
    """
    cp(adata=adata, thresh=thresh)
    sns.scatterplot(data=adata.obs, x="S_score", y="G2M_score", hue="phase", ax=ax, s=s)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    ax.vlines(thresh, ylim[0], thresh, color=color, lw=lw)
    ax.hlines(thresh, xlim[0], thresh, color=color, lw=lw)
    ax.plot(
        [thresh, max(xlim[1], ylim[1])],
        [thresh, max(xlim[1], ylim[1])],
        color=color,
        lw=lw,
    )
    ax.set_xlim(*xlim), ax.set_ylim(*ylim)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Phase", frameon=False)
    ax.set_aspect("equal")
