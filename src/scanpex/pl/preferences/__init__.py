from typing import Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Global configuration for font and path simplification
mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "path.simplify": True})


# Standard arguments for saving figures
kwarg_savefig = {
    "facecolor": "white",
    "dpi": 600,
    "bbox_inches": "tight",
    "pad_inches": 0.05,
    "transparent": True,
}

# Standard arguments for saving transparent figures
kwarg_save_transparent_fig = {
    "dpi": 600,
    "bbox_inches": "tight",
    "pad_inches": 0.05,
    "transparent": True,
}


def level_cmap(
    color: tuple = tuple(
        (3 * np.array(plt.cm.Purples(0.6)) + 2.5 * np.array([0, 0, 1, 1])) / 5
    ),
    grayscale: float = 0.2,
    gamma: float = 2,
    name: str = None,
):
    """
    Create a custom linear segmented colormap transitioning from grey to a color.

    This function generates a colormap that starts at a specific shade of grey
    and interpolates linearly to the target `color`.

    Parameters
    ----------
    color : tuple of float, optional
        The target RGBA color to transition to.
        Defaults to a specific purple-blue blend.
    grayscale : float, optional
        The intensity of the starting grey color (passed to `plt.cm.Greys`).
        0.0 is black, 1.0 is white. By default 0.2.
    gamma : float, optional
        Gamma correction factor for the colormap. By default 2.
    name : str, optional
        The name of the colormap. If None, generates a name based on the color.
        By default None.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        The generated colormap object.
    """
    cdict = {
        "red": [
            (x, y0, y1)
            for x, y0, y1 in zip(
                np.linspace(0, 1, 256),
                np.linspace(plt.cm.Greys(grayscale)[0], color[0], 256),
                np.linspace(plt.cm.Greys(grayscale)[0], color[0], 256),
            )
        ],
        "green": [
            (x, y0, y1)
            for x, y0, y1 in zip(
                np.linspace(0, 1, 256),
                np.linspace(plt.cm.Greys(grayscale)[1], color[1], 256),
                np.linspace(plt.cm.Greys(grayscale)[1], color[1], 256),
            )
        ],
        "blue": [
            (x, y0, y1)
            for x, y0, y1 in zip(
                np.linspace(0, 1, 256),
                np.linspace(plt.cm.Greys(grayscale)[2], color[2], 256),
                np.linspace(plt.cm.Greys(grayscale)[2], color[2], 256),
            )
        ],
    }

    name = name if name is not None else f"level_cmap_{color}"

    return LinearSegmentedColormap(name, cdict, gamma=gamma)


# Create a default instance for Seurat-like visualizations
seurat = level_cmap(name="seurat")


def rgba2gray(r: float, g: float, b: float, a: float) -> Tuple[float]:
    """
    Calculate the luminance of an RGBA color weighted by its alpha channel.

    Uses the standard luminance formula (0.299R + 0.587G + 0.114B) adapted
    to the provided weights.

    Parameters
    ----------
    r, g, b : float
        Red, Green, and Blue components (0.0 to 1.0).
    a : float
        Alpha component (transparency).

    Returns
    -------
    float
        The weighted luminance value.
    """
    return a * (0.299 * r + 0.578 * g + 0.114 * b)


def textcolor(c, thresh=0.4, dark=".2", light="1") -> str:
    """
    Determine whether text should be light or dark based on background brightness.

    Useful for ensuring text readability on varying background colors (e.g., in
    heatmaps or bar plots).

    Parameters
    ----------
    c : str or tuple
        The background color. Can be a string representation of a float (e.g., "0.5"
        for grayscale) or an RGBA tuple.
    thresh : float, optional
        The luminance threshold. If the calculated luminance of `c` is greater
        than this value, the background is considered "light". By default 0.4.
    dark : str, optional
        The color to return if the background is light. By default ".2" (dark grey).
    light : str, optional
        The color to return if the background is dark. By default "1" (white).

    Returns
    -------
    str
        The color code for the text (`dark` or `light`).
    """
    if isinstance(c, str):
        return dark if float(c) > thresh else light
    else:
        return dark if rgba2gray(*c) > thresh else light
