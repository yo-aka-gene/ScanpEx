from typing import Optional

import anndata as ad


def seacells(
    adata: ad.AnnData,
    seacell_size: Optional[int] = 50,
    n_SEACells: Optional[int] = None,
    n_waypoint_eigs: int = 10,
    min_iter: int = 10,
    max_iter: int = 50,
    seacell_key: str = "SEACells",
) -> None:
    """
    Compute SEACells (metacells) aggregation.

    This function is a wrapper for the SEACells algorithm.
    It identifies archetypes (metacells) based on the PCA space (`X_pca`).

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix. Must contain "X_pca" in `obsm`.
    seacell_size : int, optional (default: 50)
        Approximate number of cells per SEACell (metacell).
        Used to calculate `n_SEACells` automatically: `n_SEACells = n_obs / seacell_size`.
        If provided, this takes precedence over `n_SEACells` argument.
    n_SEACells : int, optional (default: None)
        The exact number of SEACells to identify.
        Ignored if `seacell_size` is provided.
    n_waypoint_eigs : int, optional (default: 10)
        Number of eigenvalues to use for waypoint initialization.
    min_iter : int, optional (default: 10)
        Minimum number of iterations for the model fitting.
    max_iter : int, optional (default: 50)
        Maximum number of iterations for the model fitting.
    seacell_key : str, optional (default: "SEACells")
        The key under which the hard assignments are stored in `adata.obs`.

    Returns
    -------
    None
        The result is stored in `adata.obs[seacell_key]`, containing the
        assigned SEACell ID for each single cell.

    Raises
    ------
    ImportError
        If the `SEACells` library is not installed.
    ValueError
        If neither `seacell_size` nor `n_SEACells` is provided.
    """
    try:
        import SEACells
    except ImportError:
        raise ImportError(
            "SEACells is not installed. Please run `poetry add SEACells`."
        )

    if seacell_size is not None:
        n_SEACells = int(adata.n_obs / seacell_size)

    if n_SEACells is None:
        raise ValueError("Either 'seacell_size' or 'n_SEACells' must be provided.")

    model = SEACells.core.SEACells(
        adata,
        build_kernel_on="X_pca",
        n_SEACells=n_SEACells,
        n_waypoint_eigs=n_waypoint_eigs,
    )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=min_iter, max_iter=max_iter)
    adata.obs[seacell_key] = model.get_hard_assignments()
