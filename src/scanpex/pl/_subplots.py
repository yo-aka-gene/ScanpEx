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

    elif (n is None):
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
