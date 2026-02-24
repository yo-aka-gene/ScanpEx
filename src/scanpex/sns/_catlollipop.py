from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def catlollipop(
    data: pd.DataFrame,
    x: Any = None,
    y: Any = None,
    hue: Any = None,
    palette: Any = None,
    ax: plt.Axes = None,
    **kwargs
):
    """
    Draw a categorical lollipop plot with dodging capabilities.

    This function creates a lollipop plot (a variation of a bar plot consisting of
    a segment and a dot) for categorical data. It automatically detects the
    orientation based on the numeric data type of `x` or `y`. It also supports
    grouping by a `hue` variable, which "dodges" (offsets) the lollipops
    to avoid overlap.

    Parameters
    ----------
    data : pd.DataFrame
        Input data structure.
    x, y : str, optional
        Variables that specify positions on the x and y axes. One should be
        categorical and the other numeric. The function determines orientation
        automatically.
    hue : str, optional
        Grouping variable that will produce points with different colors.
        Also used to offset (dodge) the stems for better visibility.
    palette : str, list, dict, or matplotlib.colors.Colormap, optional
        Method for choosing the colors to use when mapping the `hue` semantic.
        If None, defaults to "tab10".
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

    if pd.api.types.is_numeric_dtype(data.loc[:, x]):
        ax.set_yticks(-np.arange(data.loc[:, y].unique().size))
        ax.set_yticklabels(data.loc[:, y].unique())
        n_category = data.loc[:, y].unique().size
        horizontal = True
    elif pd.api.types.is_numeric_dtype(data.loc[:, y]):
        ax.set_xticks(-np.arange(data.loc[:, x].unique().size))
        ax.set_xticklabels(data.loc[:, x].unique())
        n_category = data.loc[:, x].unique().size
        horizontal = False
    n_marker = data.loc[:, hue].unique().size
    palette = sns.color_palette("tab10", n_marker) if palette is None else palette
    data = data.assign(__id__=list(-np.arange(n_category)) * n_marker).sort_values(
        ["__id__", hue]
    )
    data = data.assign(
        __loc__=data.__id__ - np.tile(np.linspace(-0.25, 0.25, n_marker), n_category)
    )

    args = dict(x=x, y="__loc__") if horizontal else dict(x="__loc__", y=y)
    sns.scatterplot(data=data, **args, hue=hue, palette=palette, **kwargs)

    lim = ax.get_xlim() if horizontal else ax.get_ylim()
    line_func = getattr(ax, "hlines" if horizontal else "vlines")
    set_lim = getattr(ax, "set_xlim" if horizontal else "set_ylim")

    for (i, loc), v in zip(
        enumerate(data.__loc__), data.loc[:, x if horizontal else y]
    ):
        line_func(loc, lim[0], v, color=palette[i % n_marker])

    set_lim(*lim)
    ax.legend(frameon=False)
    ax.set_ylabel("")
    return ax
