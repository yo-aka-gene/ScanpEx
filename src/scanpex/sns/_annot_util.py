def float2str(z: float, digit: int = 3) -> str:
    """
    Convert a float to a string with a specified number of decimal places.

    Parameters
    ----------
    z : floata
        The float to convert.
    digit : int, default 3
        The number of decimal places to include.

    Returns
    -------
    str
        The float as a string with the specified number of decimal places.
    """
    s = str(round(z, digit))
    subzero_digit = len(s.split(".")[-1])
    return (
        s
        if subzero_digit == digit
        else s + "".join(["0" for _ in range(digit - subzero_digit)])
    )
