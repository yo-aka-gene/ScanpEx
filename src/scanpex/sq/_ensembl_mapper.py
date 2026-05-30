from typing import List

import pandas as pd


def ensembl_mapper(
    ensembl_list: List[str],
    species: str = "human",
) -> pd.DataFrame:
    """Map a list of Ensembl gene IDs to their corresponding gene symbols.

    This function queries the MyGene.info API to translate Ensembl IDs into
    gene symbols. If a gene symbol cannot be found for a given ID, the original
    Ensembl ID is retained in the 'symbol' column. Duplicate queries are dropped,
    keeping only the first matched record.

    Parameters
    ----------
    ensembl_list : List[str]
        A list of Ensembl gene IDs to be mapped (e.g., ['ENSG00000139618']).
    species : str, optional
        The species name or taxonomy ID to restrict the query. Common values
        include 'human' or 'mouse'. Default is 'human'.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing two columns:
        - 'query': The original Ensembl ID.
        - 'symbol': The mapped gene symbol (or the original ID if unmapped).

    Raises
    ------
    ImportError
        If the `mygene` package is not installed in the current environment.

    Examples
    --------
    >>> ensembl_ids = ["ENSG00000139618", "ENSG00000157764", "INVALID_ID"]
    >>> df = ensembl_mapper(ensembl_ids, species="human")
    >>> print(df)
                 query      symbol
    0  ENSG00000139618       BRCA2
    1  ENSG00000157764        BRAF
    2       INVALID_ID  INVALID_ID
    """
    try:
        import mygene
    except ImportError:
        raise ImportError(
            "mygene is not installed. Please install it using `pip install mygene`."
        )

    mg = mygene.MyGeneInfo()
    gene_symbols = mg.querymany(
        ensembl_list,
        scopes="ensembl.gene",
        fields="symbol",
        species=species,
        as_dataframe=True,
    ).reset_index()[["query", "symbol"]]
    gene_symbols["symbol"] = gene_symbols["symbol"].fillna(gene_symbols["query"])
    return gene_symbols.drop_duplicates(subset=["query"], keep="first")
