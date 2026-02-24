from typing import Any, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def vinswarm(
    data: pd.DataFrame,
    x: Any = None,
    y: Any = None,
    hue: Any = None,
    inner: Any = None,
    cut: Union[int, float] = 0,
    alpha: Union[int, float] = 0.5,
    s: Union[int, float] = 1,
    jitter: Union[int, float] = 0.4,
    zorder: int = 0,
    ax: plt.Axes = None,
    **kwargs
) -> None:
    """
    Overlay a strip plot on a violin plot to visualize distribution and raw points.

    This function combines `seaborn.violinplot` (showing the kernel density estimate)
    and `seaborn.stripplot` (showing individual data points). It allows for
    simultaneous visualization of the underlying distribution shape and the
    sample density/outliers.

    Parameters
    ----------
    data : pd.DataFrame
        Input data structure.
    x, y : str, optional
        Variables that specify positions on the x and y axes.
    hue : str, optional
        Grouping variable that will produce elements with different colors.
    inner : str or None, optional
        Representation of the datapoints in the violin interior.
        Defaults to None (to avoid cluttering with the strip plot).
    cut : float, default 0
        Distance to extend the density past the extreme datapoints.
        Set to 0 to limit the violin range strictly within the data range.
    alpha : float, default 0.5
        Proportional opacity of the violin plot fill.
    s : float, default 1
        Size of the markers in the strip plot.
    jitter : float or bool, default 0.4
        Amount of jitter (only along the categorical axis) to apply to the
        strip plot points.
    zorder : int, default 0
        The z-order for the strip plot points. Higher values draw points
        on top of the violin.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. Otherwise, a new figure and axes
        are created.
    **kwargs
        Other keyword arguments are passed to *both* `seaborn.violinplot`
        and `seaborn.stripplot`.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes object with the plot drawn onto it.
    """
    if ax is None:
        _, ax = plt.subplots()

    sns.violinplot(
        data=data, x=x, y=y, hue=hue, inner=inner, cut=cut, alpha=alpha, **kwargs
    )
    sns.stripplot(
        data=data, x=x, y=y, hue=hue, s=s, jitter=jitter, zorder=zorder, **kwargs
    )
    return ax
