from typing import Any, List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def lollipop(
    data: pd.DataFrame,
    x: Any = None,
    y: Any = None,
    hue: Any = None,
    lim: List[float] = None,
    palette: Any = None,
    ax: plt.Axes = None,
    **kwargs
):
    """
    Draw a lollipop plot with stems and markers, supporting hue and custom limits.

    This function plots a stem (line) for each data row, extending from a baseline
    to the data point. It automatically detects the orientation based on which
    axis is numeric. It integrates with `seaborn.scatterplot` for marker
    rendering, allowing for `hue` grouping and other aesthetic customizations.

    Parameters
    ----------
    data : pd.DataFrame
        Input data structure.
    x : str, optional
        Column name for the x-axis. If numeric, a horizontal lollipop plot is
        drawn using the row index as the y-position.
    y : str, optional
        Column name for the y-axis. If numeric (and x is not), a vertical
        lollipop plot is drawn using the row index as the x-position.
    hue : str, optional
        Grouping variable passed to `seaborn.scatterplot`.
        Note: Stems are colored individually by row index using `palette`.
        If `hue` is used, ensure `palette` is compatible or accept that stems
        and markers may follow different coloring logics.
    lim : list of float, optional
        Custom limits [min, max] for the numeric axis.
        The minimum value (lim[0]) is used as the baseline for the stems.
        If None, the current axis limits are used.
    palette : list or seaborn palette, optional
        Colors to use for the stems. If None, generates a 'husl' palette
        with a distinct color for each row.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. Otherwise, a new figure and axes
        are created.
    **kwargs
        Other keyword arguments are passed through to `seaborn.scatterplot`.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes object with the plot drawn onto it.
    """
    if ax is None:
        _, ax = plt.subplots()

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    if pd.api.types.is_numeric_dtype(data.loc[:, x]):
        xlim = lim if lim is not None else xlim
        palette = (
            sns.color_palette("husl", len(data.loc[:, x]))
            if palette is None
            else palette
        )
        for i, v in enumerate(data.loc[:, x]):
            ax.hlines(i, xlim[0], v, color=palette[i])
            ax.set_xlim(*xlim) if lim is not None else None
    elif pd.api.types.is_numeric_dtype(data.loc[:, y]):
        ylim = lim if lim is not None else ylim
        palette = (
            sns.color_palette("husl", len(data.loc[:, y]))
            if palette is None
            else palette
        )
        for i, v in enumerate(data.loc[:, y]):
            ax.vlines(i, ylim[0], v, color=palette[i])
            ax.set_ylim(*ylim) if lim is not None else None

    sns.scatterplot(
        data=data, x=x, y=y, ax=ax, hue=hue, linewidth=0, palette=palette, **kwargs
    )

    ax.legend().remove()
    return ax
