from typing import Dict, Union


def multiclass_args(
    num_class: int,
    objective="multiclass",
    metric="multi_logloss",
    verbosity=-1,
    deterministic=True,
    random_seed=0,
    num_boost_round=100,
    force_col_wise=True,
) -> Dict[str, Union[str, int, bool]]:
    """
    Generate a hyperparameter dictionary for LightGBM multi-class classification.

    This helper function constructs the configuration dictionary required to
    train a LightGBM model. It sets standard defaults for reproducibility and
    logging.

    Parameters
    ----------
    num_class : int
        The number of target classes. This is a required parameter for
        multi-class objectives.
    objective : str, optional
        The learning objective. Common values include "multiclass" and
        "multiclassova". By default "multiclass".
    metric : str, optional
        The metric to be evaluated on the evaluation set. Common values
        include "multi_logloss" or "multi_error". By default "multi_logloss".
    verbosity : int, optional
        Controls the level of LightGBM logging. < 0 for fatal only, 0 for error,
        1 for info. By default -1 (silent).
    deterministic : bool, optional
        If True, ensures reproducible results. By default True.
    random_seed : int, optional
        The seed for the random number generator. By default 0.
    num_boost_round : int, optional
        The number of boosting iterations (trees) to build. By default 100.
    force_col_wise : bool, optional
        If True, forces column-wise histogram building, which can reduce memory
        usage and is generally faster on CPUs with many cores. By default True.

    Returns
    -------
    dict
        A dictionary containing the parameter keys and values ready to be
        passed to LightGBM training functions.
    """
    return dict(
        objective=objective,
        metric=metric,
        num_class=num_class,
        verbosity=verbosity,
        deterministic=deterministic,
        random_seed=random_seed,
        num_boost_round=num_boost_round,
        force_col_wise=force_col_wise,
    )


def regression_args(
    objective="regression",
    metric="l2",
    verbosity=-1,
    deterministic=True,
    random_seed=0,
    num_boost_round=100,
    force_col_wise=True,
) -> Dict[str, Union[str, int, bool]]:
    """
    Generate a hyperparameter dictionary for LightGBM regression tasks.

    Parameters
    ----------
    objective : str, optional
        The learning objective. Common values include "regression" (L2),
        "regression_l1" (L1), or "huber". By default "regression".
    metric : str, optional
        The metric to be evaluated on the evaluation set. Common values
        include "l2" (MSE), "l1" (MAE), or "rmse". By default "l2".
    verbosity : int, optional
        Controls the level of LightGBM logging. By default -1 (silent).
    deterministic : bool, optional
        If True, ensures reproducible results. By default True.
    random_seed : int, optional
        The seed for the random number generator. By default 0.
    num_boost_round : int, optional
        The number of boosting iterations (trees) to build. By default 100.
    force_col_wise : bool, optional
        If True, forces column-wise histogram building. By default True.

    Returns
    -------
    dict
        A dictionary containing the parameter keys and values ready to be
        passed to LightGBM training functions.
    """
    return dict(
        objective=objective,
        metric=metric,
        verbosity=verbosity,
        deterministic=deterministic,
        random_seed=random_seed,
        num_boost_round=num_boost_round,
        force_col_wise=force_col_wise,
    )
