import warnings

from ._violinstrip import violinstrip


def vinswarm(*args, **kwargs):
    warnings.warn(
        "`scx.sns.vinswarm` is deprecated and will be removed in a future release. "
        "Use `scx.sns.violinstrip` instead.",
        FutureWarning,
        stacklevel=2,
    )
    return violinstrip(*args, **kwargs)
