import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve


class MulticlassROC:
    """
    Compute ROC curve components and AUC scores for multi-class models.

    This class iterates through each class in the target data (One-vs-Rest approach),
    calculating the False Positive Rate (FPR), True Positive Rate (TPR), thresholds,
    and Area Under the Curve (AUC) for each.

    Attributes
    ----------
    x : list of np.ndarray
        A list containing the FPR arrays for each class.
    y : list of np.ndarray
        A list containing the TPR arrays for each class.
    thresh : list of list of float
        A list containing the thresholds used to compute the ROC curve for each class.
    auc : list of float
        A list containing the AUROC score for each class.
    """

    def __init__(self, model, x, y):
        """
        Initialize the MulticlassROC calculator.

        Parameters
        ----------
        model : object
            The trained model. It must implement a `predict` method that returns
            prediction scores or probabilities of shape (n_samples, n_classes).
        x : np.ndarray
            The feature matrix for prediction.
        y : np.ndarray
            The ground truth labels, expected to be in a multi-label format
            (one-hot encoded) of shape (n_samples, n_classes).
        """
        self.x = [
            roc_curve(y[:, i], model.predict(x)[:, i])[0] for i in range(y.shape[1])
        ]
        self.y = [
            roc_curve(y[:, i], model.predict(x)[:, i])[1] for i in range(y.shape[1])
        ]
        self.thresh = [
            roc_curve(y[:, i], model.predict(x)[:, i])[2].tolist()
            for i in range(y.shape[1])
        ]
        self.auc = [
            roc_auc_score(y[:, i], model.predict(x)[:, i], multi_class="ovr")
            for i in range(y.shape[1])
        ]


def multiclass_roc(model, x, y, ax, cmap, label_dict, minimalist: bool = False):
    """
    Plot One-vs-Rest ROC curves and the macro-average curve.

    This function visualizes the performance of a multi-class classifier by plotting
    individual ROC curves for each class and an interpolated macro-average curve.
    It also displays the baseline (random guess) and optionally the ideal curve.

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
        - If str: The name of a matplotlib colormap (e.g., "viridis").
        - If list: A list of specific colors to assign to each class.
    label_dict : list of str
        The names of the classes corresponding to the columns of `y`.
        Used for the legend.
    minimalist : bool, optional
        If True, suppresses the "ideal" curve (perfect classifier lines)
        to reduce chart clutter. By default False.

    Returns
    -------
    None
        The plot is drawn directly onto the provided `ax` object.
    """
    roc = MulticlassROC(model=model, x=x, y=y)

    for (i, x_), y_, auc, label in zip(enumerate(roc.x), roc.y, roc.auc, label_dict):
        ax.plot(
            x_,
            y_,
            label=f"{label} (AUROC:{auc.round(3)})",
            c=(
                eval(f"plt.cm.{cmap}")(i / (len(roc.auc) - 1))
                if isinstance(cmap, str)
                else cmap[i]
            ),
        )

    ax.plot(
        [0] + np.linspace(0, 1, 100).tolist(),
        [0]
        + np.stack(
            [np.interp(np.linspace(0, 1, 100), x, y) for x, y in zip(roc.x, roc.y)]
        )
        .mean(axis=0)
        .tolist(),
        c=".2",
        label=f"macro (AUROC:{np.array(roc.auc).mean().round(2)})",
    )
    baseline_name = None
    if not minimalist:
        ax.plot(
            [0, 0, 1], [0, 1, 1], linestyle=(0, (1, 2)), c=".2", label="ideal", zorder=0
        )
        baseline_name = "baseline"
    ax.plot(
        [0, 1],
        [0, 1],
        c="gray",
        label=baseline_name,
        linestyle=(0, (1, 2)),
        zorder=1,
        alpha=0.5,
    )
    ax.set(
        xlabel="false positive rate",
        ylabel="true positive rate",
        title="ROC curve (OvR)",
    )
    ax.legend(fontsize="small", frameon=False)
