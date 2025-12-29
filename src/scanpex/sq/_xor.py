from typing import Any, List, Tuple


def xor(list_a: List[Any], list_b: List[Any]) -> Tuple[List[Any], List[Any]]:
    """
    Identify elements exclusive to each of the two input lists.

    This function calculates the set difference in both directions:
    (A - B) and (B - A).

    Note
    ----
    Since this function converts inputs to sets internally:
    1. Duplicate elements in the inputs will be removed in the output.
    2. The order of elements in the output is not guaranteed.

    Parameters
    ----------
    list_a : List[Any]
        The first list to compare.
    list_b : List[Any]
        The second list to compare.

    Returns
    -------
    Tuple[List[Any], List[Any]]
        A tuple containing two lists:
        1. Elements present in `list_a` but not in `list_b`.
        2. Elements present in `list_b` but not in `list_a`.
    """
    overlap = set(list_a) & set(list_b)
    return list(set(list_a) - overlap), list(set(list_b) - overlap)
