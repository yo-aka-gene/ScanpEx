import anndata as ad
import numpy as np


def get_quantiles(data: ad.AnnData, metrics: str, by: float = 0.1, area: list = None):
    """
    Calculate quantiles for a specified metric in the observation metadata.

    Parameters
    ----------
    data : ad.AnnData
        The annotated data matrix.
    metrics : str
        The column name in `data.obs` for which to calculate quantiles.
    by : float, optional
        The step size for generating quantile thresholds if `area` is None.
        Creates a range from 0 to 1 with this step, excluding the 0th percentile.
        By default 0.1.
    area : list of float, optional
        A specific list of quantile thresholds (between 0 and 1) to compute.
        If provided, `by` is ignored. By default None.

    Returns
    -------
    list of float
        The calculated quantile values corresponding to the requested thresholds.
    """
    area = np.arange(0, 1, by)[1:] if area is None else area
    return [data.obs[metrics].quantile(v) for v in area]
