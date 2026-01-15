from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np


def subplots(
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    n: Optional[int] = None,
    base: Optional[Tuple[float, float]] = (4, 4),
    figsize: Optional[Tuple[float, float]] = None,
    vertical: bool = True,
    flatten: bool = True,
    **kwargs
) -> Tuple[plt.Figure, Union[plt.Axes, np.ndarray, List[plt.Axes]]]:
    """
    Create a figure and a set of subplots with automatic layout calculation.

    This wrapper around `matplotlib.pyplot.subplots` determines the optimal
    grid dimensions (`nrows` x `ncols`) based on the total number of plots (`n`).
    It also automatically scales the figure size based on a `base` size per
    subplot.

    Parameters
    ----------
    nrows : int, optional
        Number of rows of the subplot grid. If None and `n` is provided,
        it may be calculated automatically.
    ncols : int, optional
        Number of columns of the subplot grid. If None and `n` is provided,
        it may be calculated automatically.
    n : int, optional
        Total number of subplots required.

        - If `nrows` and `ncols` are None, the grid is calculated to be nearly
          square (controlled by `vertical`).
        - If only one of `nrows` or `ncols` is provided, the other is calculated
          to accommodate `n` plots.

    base : tuple of float, optional
        The width and height of a single subplot (width, height).
        Used to calculate `figsize` if it is not explicitly provided.
        By default (4, 4).
    figsize : tuple of float, optional
        The total figure size (width, height). If None, it is calculated as
        `(base[0] * ncols, base[1] * nrows)`.
    vertical : bool, optional
        Controls the aspect ratio logic when calculating grid dimensions from `n`.
        If True, prioritizes more rows (taller figure). If False, prioritizes
        more columns (wider figure). By default True.
    flatten : bool, optional
        If True, returns the axes as a 1D NumPy array (even if there is only one
        subplot). This ensures consistent iteration (e.g., `for ax in axs:`).
        By default True.
    **kwargs
        Additional keyword arguments passed to `matplotlib.pyplot.subplots`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure object.
    axs : matplotlib.axes.Axes or numpy.ndarray
        The axes of the subplots. If `flatten` is True, this is always a
        1D numpy array of Axes objects. Otherwise, it follows the standard
        matplotlib behavior (Axes object or array of Axes).
    """
    if (n is None) and (nrows is None) and (ncols is None):
        nrows, ncols = 1, 1

    elif isinstance(n, int):
        if (nrows is None) and (ncols is None):
            root = np.sqrt(n)
            d1 = int(np.ceil(root))
            d2 = int(np.ceil(n / d1))
            nrows = d1 if vertical else d2
            ncols = d2 if vertical else d1

        elif isinstance(nrows, int) and (ncols is None):
            ncols = int(np.ceil(n / nrows))

        elif (nrows is None) and isinstance(ncols, int):
            nrows = int(np.ceil(n / ncols))

    elif n is None:
        if isinstance(nrows, int) and (ncols is None):
            ncols = 1
        elif (nrows is None) and isinstance(ncols, int):
            nrows = 1

    if (figsize is None) and isinstance(base, (tuple, list)):
        figsize = (base[0] * ncols, base[1] * nrows)

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, **kwargs)

    if flatten:
        if isinstance(axs, np.ndarray):
            axs = axs.flatten()
        elif not isinstance(axs, list):
            axs = np.array([axs])

    return fig, axs
