from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.container import BarContainer

from ._annot_util import float2str


def barplot(
    data: pd.DataFrame,
    x: Any = None,
    y: Any = None,
    hue: Any = None,
    *,
    palette: str | list[str] | dict[Any, Any] | None = None,
    alpha: float = 0.5,
    ax: Axes | None = None,
    bar_kws: dict[str, Any] | None = None,
    annot_kws: dict[str, Any] | None = None,
) -> Axes:
    """
    Draw a bar plot with numerical annotations inside each bar.

    This function extends `seaborn.barplot` by annotating the aggregated
    value directly inside each bar. Annotations are derived from the
    rendered bars, ensuring consistency with the selected estimator,
    hue grouping, and dodge positions.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data structure.
    x, y : str, optional
        Variables defining the bar plot. By default, `x` represents
        numerical values and `y` represents categories.
    hue : str, optional
        Grouping variable used to distinguish bars by color.
    palette : str, list, or dict, optional
        Color palette passed to `seaborn.barplot`.
    alpha : float, default 0.5
        Opacity of the bars.
    ax : matplotlib.axes.Axes, optional
        Axes on which to draw the plot.
    bar_kws : dict, optional
        Additional arguments passed to `seaborn.barplot`.
    annot_kws : dict, optional
        Additional arguments passed to `Axes.bar_label`.
        The `digit` key controls decimal precision (default 3).
        By default, labels are placed 5 points inside the bar endpoints
        and displayed in white.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes containing the annotated bar plot.
    """
    if ax is None:
        _, ax = plt.subplots()

    bar_kws = {} if bar_kws is None else dict(bar_kws)
    annot_kws = {} if annot_kws is None else dict(annot_kws)

    bar_kws.setdefault("orient", "h")

    digit = annot_kws.pop("digit", 3)
    annot_kws.setdefault("label_type", "edge")
    annot_kws.setdefault("padding", -5)
    annot_kws.setdefault("color", "white")
    annot_kws.setdefault("fontsize", 10)

    n_containers = len(ax.containers)

    sns.barplot(
        data=data,
        x=x,
        y=y,
        hue=hue,
        palette=palette,
        alpha=alpha,
        ax=ax,
        **bar_kws,
    )

    for container in ax.containers[n_containers:]:
        if isinstance(container, BarContainer):
            ax.bar_label(
                container,
                fmt=lambda val: float2str(val, digit),
                **annot_kws,
            )

    return ax
