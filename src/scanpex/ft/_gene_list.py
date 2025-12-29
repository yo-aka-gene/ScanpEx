from typing import List, Dict, Tuple, Optional
import anndata as ad
import pandas as pd
from sklearn.feature_selection import RFE
from sklearn.linear_model import Lasso, LassoCV

from scanpex.sq import GeneCacheManager, gene_query
from scanpex.tl import prob_genes


class GeneList:
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
        **kwargs
    ):
        assert (database is not None) or (gene_names is not None), \
            f"Assign at least either database or gene_names"
        self.name2idx = {n: i for n, i in zip(adata.var[source_key], adata.var.index)}
        self.idx2name = {i: n for n, i in self.name2idx.items()}
        gene_names = gene_names if gene_names is not None else database[category]
        loader = GeneCacheManager()
        self.genes = loader.load(
            key=key, 
            func=(lambda gene_names: gene_names) if preset else gene_query, 
            gene_names=gene_names, 
            source=adata.var[source_key]
        )
        self.ids = [self.name2idx[g] for g in self.genes]
        self.category = category
        self.caption = key.capitalize() if caption is None else caption
        self.score_name = f"{self.caption.replace(' ', '_')}_score"
        self.score_prob_name = self.score_name + "_prob"
        
        _data = adata.copy()

        prob_genes(
            _data, 
            gene_list=self.ids, 
            score_name=self.score_name, 
            copy=False,
            **kwargs
        )

        self.data = _data[:, self.ids].copy()


    def _get_aggregated_df(
        self,
        group_key: Optional[str],
        score_name: str
    ) -> pd.DataFrame:
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
        **kwargs
    ) -> Tuple[List[str], List[str]]:
        score_name = self.score_name if use_raw_score else self.score_prob_name
        data = self._get_aggregated_df(group_key, score_name)
        self.selected_ids = data.corr(**kwargs)["score"].drop("score").sort_values(
            ascending=False
        ).iloc[:n_top].index.tolist()
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
        **kwargs
    ) -> Tuple[List[str], List[str]]:
        score_name = self.score_name if use_raw_score else self.score_prob_name
        data = self._get_aggregated_df(group_key, score_name)
        X = data.loc[:, self.ids]
        y = data.loc[:, "score"].values.ravel()
        lasso_cv = LassoCV(cv=n_cv, random_state=random_state).fit(X, y)
        best_alpha = lasso_cv.alpha_
        selector = RFE(
            estimator=Lasso(alpha=best_alpha),
            n_features_to_select=n_top,
            step=step
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
        use_gene_name: bool = True
    ) -> pd.DataFrame:
        if not hasattr(self, 'selected_ids'):
             raise RuntimeError("Run select_correlated_genes or select_independent_genes first.")

        score_name = self.score_name if use_raw_score else self.score_prob_name
        ids = self.selected_ids + ["score"] if with_score else self.selected_ids
        genes = self.selected_genes + [score_name] if with_score else self.selected_genes
        data = self._get_aggregated_df(group_key, score_name).loc[:, ids]

        if use_gene_name:
            data.columns = genes

        return data
