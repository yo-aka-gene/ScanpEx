from typing import Any

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import pandas as pd
import seaborn as sns


def lollipop(
    data: pd.DataFrame,
    x: Any = None,
    y: Any = None,
    hue: Any = None,
    lim: list[float] | None = None,
    palette: Any = None,
    ax: Axes | None = None,
    *,
    size: Any = None,
    style: Any = None,
    sizes: Any = None,
    markers: Any = None,
    legend: str | bool = "auto",
    marker_kws: dict[str, Any] | None = None,
    stem_kws: dict[str, Any] | None = None,
    **kwargs,
) -> Axes:
    """
    Draw a lollipop plot with configurable marker aesthetics.

    This function draws a stem for each observation from a common baseline to
    its corresponding value and overlays the stem endpoints using
    `seaborn.scatterplot`. Marker color, size, and style can therefore encode
    additional variables through the standard Seaborn semantic mappings.

    The plot orientation is determined automatically from the numeric axis. If
    `x` is numeric, a horizontal lollipop plot is drawn. Otherwise, if `y` is
    numeric, a vertical lollipop plot is drawn.

    Parameters
    ----------
    data : pd.DataFrame
        Input data structure.
    x, y : str, optional
        Variables that specify positions on the x and y axes. One axis should
        contain numeric values defining the stem lengths.
    hue : str, optional
        Variable mapped to marker color. When specified, the corresponding
        stems use the same hue-color mapping as the endpoint markers.
    lim : list of float, optional
        Limits for the numeric axis. The lower limit is used as the common
        baseline from which the stems are drawn. If None, the current lower
        axis limit is used.
    palette : str, sequence, dict, or seaborn palette, optional
        Color palette used for the plot. When `hue` is specified, the palette
        defines the mapping from hue levels to colors and is shared by the
        markers and stems. Without `hue`, the palette is used to color stems.
    ax : matplotlib.axes.Axes, optional
        Pre-existing Axes on which to draw the plot. If None, a new figure and
        Axes are created.
    size : str, optional
        Variable mapped to marker size through `seaborn.scatterplot`.
    style : str, optional
        Variable mapped to marker style through `seaborn.scatterplot`.
    sizes : list, tuple, dict, or tuple of float, optional
        Specification controlling the marker-size mapping passed to
        `seaborn.scatterplot`.
    markers : bool, list, or dict, optional
        Specification controlling the marker-style mapping passed to
        `seaborn.scatterplot`.
    legend : {"auto", "brief", "full"} or bool, default "auto"
        Controls generation of the semantic legend by `seaborn.scatterplot`.
        Set to False to suppress the legend.
    marker_kws : dict, optional
        Additional keyword arguments passed only to `seaborn.scatterplot`.
    stem_kws : dict, optional
        Additional keyword arguments passed only to the stem layer drawn with
        `Axes.hlines` or `Axes.vlines`.
    **kwargs
        Additional keyword arguments passed to `seaborn.scatterplot`.

        This argument is retained for backward compatibility. New code should
        prefer `marker_kws` for marker-specific customization.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes object containing the lollipop plot.

    Raises
    ------
    ValueError
        If neither `x` nor `y` refers to a numeric column.

    Notes
    -----
    Explicit semantic arguments such as `hue`, `size`, and `style` take
    precedence over values supplied through `marker_kws` or legacy `**kwargs`.

    When `hue` is specified, a common categorical color mapping is constructed
    and used for both the markers and their corresponding stems.
    """
    if ax is None:
        _, ax = plt.subplots()

    marker_kws = {} if marker_kws is None else dict(marker_kws)
    stem_kws = {} if stem_kws is None else dict(stem_kws)

    # Preserve legacy **kwargs while allowing marker_kws to override them.
    scatter_kws = dict(kwargs)
    scatter_kws.update(marker_kws)
    scatter_kws.setdefault("linewidth", 0)

    x_is_numeric = (
        x is not None
        and pd.api.types.is_numeric_dtype(data.loc[:, x])
    )
    y_is_numeric = (
        y is not None
        and pd.api.types.is_numeric_dtype(data.loc[:, y])
    )

    if x_is_numeric:
        horizontal = True
        values = data.loc[:, x]
        axis_lim = lim if lim is not None else ax.get_xlim()
    elif y_is_numeric:
        horizontal = False
        values = data.loc[:, y]
        axis_lim = lim if lim is not None else ax.get_ylim()
    else:
        raise ValueError(
            "Either `x` or `y` must refer to a numeric column."
        )

    # Construct colors for the stem layer.
    scatter_palette = palette

    if hue is not None:
        hue_levels = list(pd.unique(data.loc[:, hue]))

        if isinstance(palette, dict):
            hue_palette = palette
        else:
            colors = sns.color_palette(
                palette,
                n_colors=len(hue_levels),
            )
            hue_palette = dict(zip(hue_levels, colors))

        stem_colors = [
            hue_palette[value]
            for value in data.loc[:, hue]
        ]
        scatter_palette = hue_palette

    else:
        if palette is None:
            stem_colors = sns.color_palette(
                "husl",
                n_colors=len(data),
            )
        else:
            stem_colors = sns.color_palette(
                palette,
                n_colors=len(data),
            )

    # Draw stems.
    for i, (value, color) in enumerate(
        zip(values, stem_colors)
    ):
        line_kws = dict(stem_kws)
        line_kws.setdefault("color", color)

        if horizontal:
            ax.hlines(
                i,
                axis_lim[0],
                value,
                **line_kws,
            )
        else:
            ax.vlines(
                i,
                axis_lim[0],
                value,
                **line_kws,
            )

    if lim is not None:
        if horizontal:
            ax.set_xlim(*lim)
        else:
            ax.set_ylim(*lim)

    # Explicit semantic arguments take precedence over generic kwargs.
    scatter_kws["hue"] = hue
    scatter_kws["size"] = size
    scatter_kws["style"] = style
    scatter_kws["legend"] = legend

    if hue is not None:
        scatter_kws["palette"] = scatter_palette
    else:
        scatter_kws.pop("palette", None)

    if size is not None:
        scatter_kws["sizes"] = sizes
    else:
        scatter_kws.pop("sizes", None)

    if style is not None:
        scatter_kws["markers"] = markers
    else:
        scatter_kws.pop("markers", None)

    sns.scatterplot(
        data=data,
        x=x,
        y=y,
        ax=ax,
        **scatter_kws,
    )

    return ax
    