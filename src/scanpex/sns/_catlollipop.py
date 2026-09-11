from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes


def catlollipop(
    data: pd.DataFrame,
    x: Any = None,
    y: Any = None,
    hue: Any = None,
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
    dodge: float = 0.25,
    **kwargs,
) -> Axes:
    """
    Draw a categorical lollipop plot with dodged marker groups.

    This function draws a stem and endpoint marker for each observation while
    arranging observations along a categorical axis. When `hue` is specified,
    groups within each category are offset along the categorical axis to avoid
    overlap.

    Marker rendering is delegated to `seaborn.scatterplot`, allowing marker
    color, size, and style to encode additional variables.

    Parameters
    ----------
    data : pd.DataFrame
        Input data structure.
    x, y : str, optional
        Variables that specify positions on the x and y axes. Exactly one
        should be numeric and the other categorical. The orientation is
        determined automatically from the numeric variable.
    hue : str, optional
        Grouping variable mapped to marker and stem color. Hue groups are also
        dodged within each categorical position.
    palette : str, sequence, dict, or seaborn palette, optional
        Color palette used for the `hue` mapping. The same color mapping is
        applied to markers and their corresponding stems.
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
    dodge : float, default 0.25
        Maximum offset from the center of each categorical position used to
        separate hue groups. Has no effect when `hue` is None.
    **kwargs
        Additional keyword arguments passed to `seaborn.scatterplot`.

        This argument is retained for backward compatibility. New code should
        prefer `marker_kws` for marker-specific customization.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes object containing the categorical lollipop plot.

    Raises
    ------
    ValueError
        If both `x` and `y` are numeric, neither is numeric, or `dodge` is
        negative.

    Notes
    -----
    Unlike implementations that assume a complete Cartesian product of
    category and hue levels, positions are calculated independently for each
    observation. Missing category-hue combinations therefore do not affect the
    positions of other observations.

    When `hue` is specified, it has two roles: it controls marker and stem
    color and also determines the dodge group within each category.
    """
    if ax is None:
        _, ax = plt.subplots()

    if dodge < 0:
        raise ValueError("`dodge` must be non-negative.")

    marker_kws = {} if marker_kws is None else dict(marker_kws)
    stem_kws = {} if stem_kws is None else dict(stem_kws)

    # Preserve legacy **kwargs while allowing marker_kws to override them.
    scatter_kws = dict(kwargs)
    scatter_kws.update(marker_kws)
    scatter_kws.setdefault("linewidth", 0)

    x_is_numeric = x is not None and pd.api.types.is_numeric_dtype(data.loc[:, x])
    y_is_numeric = y is not None and pd.api.types.is_numeric_dtype(data.loc[:, y])

    if x_is_numeric and not y_is_numeric:
        horizontal = True
        numeric = x
        category = y
    elif y_is_numeric and not x_is_numeric:
        horizontal = False
        numeric = y
        category = x
    else:
        raise ValueError("Exactly one of `x` and `y` must refer to a numeric column.")

    plot_data = data.copy()

    # Map each category to a fixed base position.
    categories = list(pd.unique(plot_data.loc[:, category]))
    category_positions = {value: -i for i, value in enumerate(categories)}

    loc_col = "__scanpex_loc__"
    while loc_col in plot_data.columns:
        loc_col = f"_{loc_col}"

    plot_data.loc[:, loc_col] = (
        plot_data.loc[:, category].map(category_positions).astype(float)
    )

    # Construct hue-specific dodge offsets and colors.
    hue_palette = None

    if hue is not None:
        hue_levels = list(pd.unique(plot_data.loc[:, hue]))
        n_hue = len(hue_levels)

        if n_hue == 1:
            offsets = np.array([0.0])
        else:
            offsets = np.linspace(-dodge, dodge, n_hue)

        hue_offsets = dict(zip(hue_levels, offsets))

        plot_data.loc[:, loc_col] += (
            plot_data.loc[:, hue].map(hue_offsets).astype(float)
        )

        if isinstance(palette, dict):
            missing = [level for level in hue_levels if level not in palette]
            if missing:
                raise ValueError(
                    "`palette` does not define colors for all hue levels: " f"{missing}"
                )

            hue_palette = {level: palette[level] for level in hue_levels}

        else:
            colors = sns.color_palette(
                palette or "tab10",
                n_colors=n_hue,
            )
            hue_palette = dict(zip(hue_levels, colors))

    # Configure scatter semantics.
    scatter_kws["hue"] = hue
    scatter_kws["size"] = size
    scatter_kws["style"] = style
    scatter_kws["legend"] = legend

    if hue is not None:
        scatter_kws["palette"] = hue_palette
    else:
        scatter_kws.pop("palette", None)

        if palette is None:
            default_color = sns.color_palette(
                "tab10",
                n_colors=1,
            )[0]
        elif isinstance(palette, dict):
            if not palette:
                raise ValueError("`palette` must not be empty.")
            default_color = next(iter(palette.values()))
        else:
            default_color = sns.color_palette(
                palette,
                n_colors=1,
            )[0]

        scatter_kws.setdefault("color", default_color)
        stem_kws.setdefault("color", default_color)

    if size is not None:
        scatter_kws["sizes"] = sizes
    else:
        scatter_kws.pop("sizes", None)

    if style is not None:
        scatter_kws["markers"] = markers
    else:
        scatter_kws.pop("markers", None)

    scatter_args = (
        {"x": numeric, "y": loc_col} if horizontal else {"x": loc_col, "y": numeric}
    )

    # Draw markers first so that the numeric axis limits are established.
    sns.scatterplot(
        data=plot_data,
        ax=ax,
        **scatter_args,
        **scatter_kws,
    )

    numeric_lim = ax.get_xlim() if horizontal else ax.get_ylim()

    # Draw stems using the same hue-color mapping as the markers.
    line_func = ax.hlines if horizontal else ax.vlines

    for _, row in plot_data.iterrows():
        line_kws = dict(stem_kws)

        if hue is not None:
            line_kws["color"] = hue_palette[row[hue]]

        line_func(
            row[loc_col],
            numeric_lim[0],
            row[numeric],
            **line_kws,
        )

    # Restore numeric limits after adding the stems.
    if horizontal:
        ax.set_xlim(*numeric_lim)
        ax.set_yticks(list(category_positions.values()))
        ax.set_yticklabels(categories)
        ax.set_ylabel(category)
    else:
        ax.set_ylim(*numeric_lim)
        ax.set_xticks(list(category_positions.values()))
        ax.set_xticklabels(categories)
        ax.set_xlabel(category)

    if ax.get_legend() is not None:
        ax.get_legend().set_frame_on(False)

    return ax
