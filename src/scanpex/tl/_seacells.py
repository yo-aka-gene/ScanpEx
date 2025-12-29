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
