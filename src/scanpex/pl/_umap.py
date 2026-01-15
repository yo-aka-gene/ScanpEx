import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

sc.set_figure_params(scanpy=False, vector_friendly=True, dpi_save=600)


def umap(adata, title: str = None, **kwargs):
    """
    Generate a UMAP plot with customized aesthetics and legend placement.

    This function wraps `scanpy.pl.umap` to provide a cleaner look by removing
    axis spines/ticks, enforcing an equal aspect ratio, and positioning the
    legend (for categorical data) or colorbar (for continuous data) outside
    the plotting area.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix containing UMAP coordinates.
    title : str, optional
        The title to display for the legend or colorbar. If None, the
        capitalized name of the `color` column is used. By default None.
    **kwargs
        Additional keyword arguments passed to `scanpy.pl.umap`.
        Notable arguments handled specifically by this wrapper:

        - ``color`` (str): Key in `adata.obs` to color observations by.
        - ``ax`` (matplotlib.axes.Axes): The axis to plot on. If not provided,
          a new figure and axis are created.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes) or None
        Returns a tuple of `(fig, ax)` if a new figure was created (i.e., `ax`
        was not provided in kwargs).
        Returns `None` if an existing `ax` was passed (modifies the axis in-place).
    """
    ax = kwargs.get("ax", None)
    create_new_ax = True if ax is None else False
    if ax is None:
        fig, ax = plt.subplots()

    color = kwargs.get("color", None)

    kwargs = {**kwargs, **dict(ax=ax, show=False, colorbar_loc=None)}

    sc.pl.umap(adata, **kwargs)

    if (color is not None) and isinstance(adata.obs[color].dtype, pd.CategoricalDtype):
        ax.legend(
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            title=color.capitalize() if title is None else title,
            frameon=False,
        )

    elif (color is not None) and not isinstance(
        adata.obs[color].dtype, pd.CategoricalDtype
    ):
        mappable = ax.collections[0]
        cax = ax.inset_axes([1.02, 0.1, 0.04, 0.3])
        cbar = plt.colorbar(mappable, cax=cax)
        cbar.set_label(
            color.capitalize() if title is None else title, rotation=90, labelpad=5
        )
        cbar.outline.set_visible(False)

    ax.axis("off")
    ax.set_aspect("equal")
    ax.set(title="")
    return (fig, ax) if create_new_ax else None
