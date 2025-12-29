import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve


class MulticlassROC:
    def __init__(self, model, x, y):
        self.x = [
            roc_curve(y[:, i], model.predict(x)[:, i])[0] for i in range(y.shape[1])
        ]
        self.y = [
            roc_curve(y[:, i], model.predict(x)[:, i])[1] for i in range(y.shape[1])
        ]
        self.thresh = [
            roc_curve(y[:, i], model.predict(x)[:, i])[2].tolist() for i in range(y.shape[1])
        ]
        self.auc = [
            roc_auc_score(y[:, i], model.predict(x)[:, i], multi_class='ovr') for i in range(y.shape[1])
        ]


def multiclass_roc(model, x, y, ax, cmap, label_dict, minimalist: bool = False):    
    roc = MulticlassROC(model=model, x=x, y=y)
    
    for (i, x_), y_, auc, label in zip(enumerate(roc.x), roc.y, roc.auc, label_dict):
        ax.plot(
            x_, y_, 
            label=f"{label} (AUROC:{auc.round(3)})", 
            c=eval(f"plt.cm.{cmap}")(i/(len(roc.auc) - 1)) if isinstance(cmap, str) else cmap[i]
        )
    
    ax.plot(
        [0] + np.linspace(0, 1, 100).tolist(), 
        [0] + np.stack([
            np.interp(np.linspace(0, 1, 100), x, y) for x, y in zip(roc.x, roc.y)
        ]).mean(axis=0).tolist(),
        c=".2",
        label=f"macro (AUROC:{np.array(roc.auc).mean().round(2)})"
    )
    baseline_name = None
    if not minimalist: 
        ax.plot([0, 0, 1], [0, 1, 1], linestyle=(0, (1, 2)), c=".2", label="ideal", zorder=0)
        baseline_name = "baseline"
    ax.plot([0, 1], [0, 1], c="gray", label=baseline_name, linestyle=(0, (1, 2)), zorder=1, alpha=0.5)
    ax.set(xlabel="false positive rate", ylabel="true positive rate", title="ROC curve (OvR)")
    ax.legend(fontsize="small", frameon=False)
