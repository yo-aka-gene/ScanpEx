from typing import List

import numpy as np
import pandas as pd


def gene_query(
    gene_names: list,
    source: list,
    species: str = "human",
    logging: bool = True,
    unique: bool = True,
    sort: bool = False,
    keep_unmapped: bool = False,
) -> List[str]:
    """
    Map gene names (symbols or aliases) to a target source list (e.g., `adata.var_names`).

    This function uses MyGene.info to resolve gene aliases. It checks if the queried
    gene or its aliases exist in the provided `source` list. If a match is found,
    the gene name as it appears in `source` is returned.

    Parameters
    ----------
    gene_names : list
        List of gene names or aliases to query.
    source : list
        The target list of valid gene names (e.g., `adata.var_names`).
        The function checks if the queried genes exist in this list.
    species : str, optional (default: "human")
        Species to query in MyGene.info (e.g., "human", "mouse").
    logging : bool, optional (default: True)
        If True, prints the number of mapped genes and missing queries.
    unique : bool, optional (default: True)
        If True, returns a sorted list of unique gene names.
        If False, allows duplicates and maintains the original query order.
    sort : bool, optional (default: False)
        If True, sorts the returned list of genes alphanumerically.
    keep_unmapped : bool, optional (default: False)
        If True, includes unmapped gene names in the returned list.
        If False, omits unmapped genes.

    Returns
    -------
    List[str]
        A list of gene names that were successfully mapped to the `source`.

    Raises
    ------
    ImportError
        If the `mygene` library is not installed.
    """
    try:
        import mygene
    except ImportError:
        raise ImportError(
            "mygene is not installed. Please install it using `pip install mygene`."
        )

    mg = mygene.MyGeneInfo()
    res = mg.querymany(
        gene_names,
        scopes="symbol,alias,ensembl.gene",
        fields="symbol,alias",
        species=species,
        as_dataframe=True,
    )

    source_set = set(source)
    mapping_dict = {}
    found_count = 0

    for query in np.unique(gene_names):
        if query not in res.index:
            mapping_dict[query] = query if keep_unmapped else None
            continue

        match_rows = res.loc[[query]]
        candidates = []

        for _, row in match_rows.iterrows():
            if not pd.isna(row.get("symbol")):
                candidates.append(row["symbol"])
            aliases = row.get("alias")
            if isinstance(aliases, list):
                candidates.extend(aliases)
            elif isinstance(aliases, str):
                candidates.append(aliases)
            candidates.append(query)

        candidates = list(set(candidates))
        match_found = False

        for cand in candidates:
            if cand in source_set:
                mapping_dict[query] = cand
                match_found = True
                found_count += 1
                break

        if not match_found:
            if logging:
                print(f"Not found in source: {query} (Candidates: {candidates})")
            mapping_dict[query] = query if keep_unmapped else None

    if logging:
        n_total = len(np.unique(gene_names))
        print(f"[{found_count}/{n_total}] queries mapped to the source.")

    if unique:
        final_genes = list(set([v for v in mapping_dict.values() if v is not None]))
        return sorted(final_genes) if sort else final_genes

    else:
        final_genes = [
            mapping_dict[query]
            for query in gene_names
            if mapping_dict[query] is not None
        ]
        return sorted(final_genes) if sort else final_genes
