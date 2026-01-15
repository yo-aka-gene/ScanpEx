import numpy as np
from sklearn.metrics import average_precision_score, precision_recall_curve


class MulticlassPR:
    """
    Compute Precision-Recall curve components and AP scores for multi-class models.

    This class calculates the Recall (x-axis), Precision (y-axis), thresholds,
    and Average Precision (AP) for each class individually. It also calculates
    the baseline precision (prevalence) for each class.

    Attributes
    ----------
    x : list of np.ndarray
        A list containing the Recall values for each class.
    y : list of np.ndarray
        A list containing the Precision values for each class.
    thresh : list of list of float
        A list containing the thresholds used to compute the PR curve.
    ap : list of float
        A list containing the Average Precision (AP) score for each class.
    base : list of float
        A list containing the prevalence (fraction of positives) for each class.
        This represents the performance of a random classifier ("no-skill" line).
    """

    def __init__(self, model, x, y):
        """
        Initialize the MulticlassPR calculator.

        Parameters
        ----------
        model : object
            The trained model. Its `predict` method should return scores or
            probabilities suitable for curve calculation.
        x : np.ndarray
            The feature matrix.
        y : np.ndarray
            The ground truth labels (one-hot encoded) of shape (n_samples, n_classes).
        """
        self.x = [
            precision_recall_curve(y[:, i], model.predict(x)[:, i])[1]
            for i in range(y.shape[1])
        ]
        self.y = [
            precision_recall_curve(y[:, i], model.predict(x)[:, i])[0]
            for i in range(y.shape[1])
        ]
        self.thresh = [
            precision_recall_curve(y[:, i], model.predict(x)[:, i])[2].tolist()
            for i in range(y.shape[1])
        ]
        self.ap = [
            average_precision_score(y[:, i], model.predict(x)[:, i])
            for i in range(y.shape[1])
        ]
        self.base = [len(y[:, i][y[:, i] == 1]) / len(y) for i in range(y.shape[1])]


def multiclass_pr(model, x, y, ax, cmap, label_dict, minimalist: bool = False):
    """
    Plot One-vs-Rest Precision-Recall curves and the micro-average curve.

    This function visualizes the trade-off between Precision and Recall.
    It plots:
    1. Individual curves for each class.
    2. A "micro-average" curve calculated by aggregating all class predictions.
    3. Baseline lines representing the prevalence of each class.
    4. An optional "ideal" line (Precision=1.0).

    Parameters
    ----------
    model : object
        The trained classification model.
    x : np.ndarray
        The input features.
    y : np.ndarray
        The ground truth labels (one-hot encoded).
    ax : matplotlib.axes.Axes
        The axis on which to draw the plot.
    cmap : str or list
        The color mapping strategy.
        - If str: Name of a matplotlib colormap.
        - If list: Specific colors for each class.
    label_dict : list of str
        The names of the classes corresponding to the columns of `y`.
    minimalist : bool, optional
        If True, hides the "ideal" line (perfect precision) to reduce clutter.
        By default False.

    Returns
    -------
    None
        The plot is drawn directly onto the provided `ax` object.
    """
    pr = MulticlassPR(model=model, x=x, y=y)

    for (i, x_), y_, ap, label, base in zip(
        enumerate(pr.x), pr.y, pr.ap, label_dict, pr.base
    ):
        ax.plot(
            x_,
            y_,
            label=f"{label} (AP:{ap.round(3)})",
            c=(
                eval(f"plt.cm.{cmap}")(i / (len(pr.ap) - 1))
                if isinstance(cmap, str)
                else cmap[i]
            ),
        )
        ax.plot(
            [0, 1], [base, base], c="gray", linestyle=(0, (1, 2)), zorder=1, alpha=0.5
        )

    p_avg, r_avg, thr_avg = precision_recall_curve(y.ravel(), model.predict(x).ravel())

    ax.plot(
        r_avg,
        p_avg,
        c=".2",
        label=f"micro (AP:{average_precision_score(y, model.predict(x)).round(2)})",
    )
    base = np.array(pr.base).mean()

    baseline_name = None
    if not minimalist:
        ax.plot([0, 1], [1, 1], linestyle=(0, (1, 2)), c=".2", label="ideal", zorder=2)
        baseline_name = "baselines"
    ax.plot(
        [0, 1],
        [base, base],
        c="gray",
        linestyle=(0, (1, 2)),
        label=baseline_name,
        zorder=1,
        alpha=0.5,
    )
    ax.set(xlabel="recall", ylabel="precision", title="PR curve (OvR)")
    ax.legend(fontsize="small", frameon=False)
