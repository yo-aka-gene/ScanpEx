from typing import Any, Tuple, List


def xor(list_a: List[Any], list_b: List[Any]) -> Tuple[List[Any], List[Any]]:
    overlap = set(list_a) & set(list_b)
    return list(set(list_a) - overlap), list(set(list_b) - overlap)
