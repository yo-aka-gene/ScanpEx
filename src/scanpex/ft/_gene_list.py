from typing import Dict, List, Optional, Tuple

import anndata as ad
import pandas as pd
from sklearn.feature_selection import RFE
from sklearn.linear_model import Lasso, LassoCV

from scanpex.sq import GeneCacheManager, gene_query
from scanpex.tl import prob_genes


class GeneList:
    """
    Manages a specific set of genes for scoring, aggregation, and feature selection.

    This class handles the retrieval of gene lists (via cache or direct input),
    calculates module scores (using `prob_genes`), and creates a subsetted
    AnnData object. It provides methods to select representative genes based on
    correlation or independence (Lasso-based).

    Attributes
    ----------
    genes : list of str
        The list of gene symbols (names).
    ids : list of str
        The list of gene IDs (indices in `adata.var`).
    data : ad.AnnData
        A subset of the original AnnData containing only the selected genes
        and calculated scores.
    score_name : str
        The key for the raw module score in `data.obs`.
    score_prob_name : str
        The key for the probabilistic (transformed) score in `data.obs`.
    category : str
        The category name used to retrieve genes from the database.
    caption : str
        The display caption for the gene list.
    """

    def __init__(
        self,
        adata: ad.AnnData,
        key: str,
        category: str,
        database: Dict[str, str] = None,
        gene_names: List[str] = None,
        preset: bool = False,
        source_key: str = "gene_name",
        caption: str = None,
        **kwargs,
    ):
        """
        Initialize the GeneList object and calculate scores.

        Parameters
        ----------
        adata : ad.AnnData
            The annotated data matrix.
        key : str
            A unique identifier for caching the gene list.
        category : str
            The key to look up in `database` if `gene_names` is not provided.
        database : dict, optional
            A dictionary mapping categories to lists of genes. Required if
            `gene_names` is None.
        gene_names : list of str, optional
            An explicit list of gene names. If provided, `database` is ignored.
        preset : bool, optional
            If True, assumes the provided names are final and skips the query
            step. By default False.
        source_key : str, optional
            The column in `adata.var` containing gene symbols.
            By default "gene_name".
        caption : str, optional
            A display name for the score. If None, derived from `key`.
            By default None.
        **kwargs
            Additional keyword arguments passed to `scanpex.tl.prob_genes`
            for score calculation.
        """
        assert (database is not None) or (
            gene_names is not None
        ), "Assign at least either database or gene_names"
        self.name2idx = {n: i for n, i in zip(adata.var[source_key], adata.var.index)}
        self.idx2name = {i: n for n, i in self.name2idx.items()}
        gene_names = gene_names if gene_names is not None else database[category]
        loader = GeneCacheManager()
        self.genes = loader.load(
            key=key,
            func=(lambda gene_names: gene_names) if preset else gene_query,
            gene_names=gene_names,
            source=adata.var[source_key],
        )
        self.ids = [self.name2idx[g] for g in self.genes]
        self.category = category
        self.caption = key.capitalize() if caption is None else caption
        self.score_name = f"{self.caption.replace(' ', '_')}_score"
        self.score_prob_name = self.score_name + "_prob"

        _data = adata.copy()

        prob_genes(
            _data, gene_list=self.ids, score_name=self.score_name, copy=False, **kwargs
        )

        self.data = _data[:, self.ids].copy()

    def _get_aggregated_df(
        self, group_key: Optional[str], score_name: str
    ) -> pd.DataFrame:
        """
        Aggregate expression data and scores by a grouping key.

        Parameters
        ----------
        group_key : str or None
            The key in `obs` to group by. If None, returns individual cell data.
        score_name : str
            The name of the score column to include.

        Returns
        -------
        pd.DataFrame
            The aggregated (mean) DataFrame.
        """
        df = self.data[:, self.ids].to_df()

        if score_name in self.data.obs:
            df["score"] = self.data.obs[score_name].values

        if group_key is not None:
            if group_key not in self.data.obs:
                raise KeyError(f"Key '{group_key}' not found in adata.obs")

            df["_group"] = self.data.obs[group_key].values

            return df.groupby("_group", observed=False).mean()

        return df

    def select_correlated_genes(
        self,
        n_top: int,
        group_key: Optional[str] = "SEACells",
        use_raw_score: bool = False,
        **kwargs,
    ) -> Tuple[List[str], List[str]]:
        """
        Select genes most highly correlated with the module score.

        This method aggregates data by `group_key` (e.g., metacells) and
        computes the correlation between each gene's expression and the
        module score.

        Parameters
        ----------
        n_top : int
            The number of top correlated genes to select.
        group_key : str, optional
            The key in `obs` used for aggregation before correlation.
            By default "SEACells".
        use_raw_score : bool, optional
            If True, uses the raw score (`score_name`).
            If False, uses the probabilistic score (`score_prob_name`).
            By default False.
        **kwargs
            Additional arguments passed to `pandas.DataFrame.corr`.

        Returns
        -------
        tuple of (list of str, list of str)
            A tuple containing (selected_gene_names, selected_gene_ids).
        """
        score_name = self.score_name if use_raw_score else self.score_prob_name
        data = self._get_aggregated_df(group_key, score_name)
        self.selected_ids = (
            data.corr(**kwargs)["score"]
            .drop("score")
            .sort_values(ascending=False)
            .iloc[:n_top]
            .index.tolist()
        )
        self.selected_genes = [self.idx2name[g] for g in self.selected_ids]
        return (self.selected_genes, self.selected_ids)

    def select_independent_genes(
        self,
        n_top: int,
        group_key: Optional[str] = "SEACells",
        use_raw_score: bool = False,
        n_cv: int = 5,
        step: float = 10,
        random_state: int = 0,
        **kwargs,
    ) -> Tuple[List[str], List[str]]:
        """
        Select a subset of genes that independently predict the module score.

        This method uses a two-step process:

        1. LassoCV to determine the optimal regularization parameter (alpha).
        2. Recursive Feature Elimination (RFE) with Lasso to select the
           top `n_top` features.

        Parameters
        ----------
        n_top : int
            The number of features to select.
        group_key : str, optional
            The key in `obs` used for aggregation. By default "SEACells".
        use_raw_score : bool, optional
            If True, uses the raw score as the target variable.
            By default False.
        n_cv : int, optional
            Number of folds for LassoCV. By default 5.
        step : float or int, optional
            Number of features to remove at each iteration of RFE. By default 10.
        random_state : int, optional
            Seed for reproducibility. By default 0.
        **kwargs
            (Currently unused, but accepted for compatibility).

        Returns
        -------
        tuple of (list of str, list of str)
            A tuple containing (selected_gene_names, selected_gene_ids).
        """
        score_name = self.score_name if use_raw_score else self.score_prob_name
        data = self._get_aggregated_df(group_key, score_name)
        X = data.loc[:, self.ids]
        y = data.loc[:, "score"].values.ravel()
        lasso_cv = LassoCV(cv=n_cv, random_state=random_state).fit(X, y)
        best_alpha = lasso_cv.alpha_
        selector = RFE(
            estimator=Lasso(alpha=best_alpha), n_features_to_select=n_top, step=step
        )
        selector.fit(X, y)
        self.selected_ids = [g for g, s in zip(self.ids, selector.support_) if s]
        self.selected_genes = [self.idx2name[g] for g in self.selected_ids]
        return (self.selected_genes, self.selected_ids)

    def get_matrix(
        self,
        group_key: Optional[str] = "SEACells",
        with_score: bool = False,
        use_raw_score: bool = False,
        use_gene_name: bool = True,
    ) -> pd.DataFrame:
        """
        Retrieve a DataFrame of the selected genes.

        Requires running `select_correlated_genes` or `select_independent_genes`
        first.

        Parameters
        ----------
        group_key : str, optional
            The key for aggregation. By default "SEACells".
        with_score : bool, optional
            If True, includes the score column in the output. By default False.
        use_raw_score : bool, optional
            If True, uses the raw score instead of the probabilistic score
            (relevant only if `with_score` is True). By default False.
        use_gene_name : bool, optional
            If True, sets columns to gene symbols. If False, uses gene IDs.
            By default True.

        Returns
        -------
        pd.DataFrame
            The data matrix of selected genes.

        Raises
        ------
        RuntimeError
            If no genes have been selected yet.
        """
        if not hasattr(self, "selected_ids"):
            raise RuntimeError(
                "Run select_correlated_genes or select_independent_genes first."
            )

        score_name = self.score_name if use_raw_score else self.score_prob_name
        ids = self.selected_ids + ["score"] if with_score else self.selected_ids
        genes = (
            self.selected_genes + [score_name] if with_score else self.selected_genes
        )
        data = self._get_aggregated_df(group_key, score_name).loc[:, ids]

        if use_gene_name:
            data.columns = genes

        return data
