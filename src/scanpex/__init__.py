from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scanpex")
except PackageNotFoundError:
    __version__ = "unknown"

from . import ft, ml, pl, pp, sq, tl

__all__ = [
    "ft",
    "ml",
    "pl",
    "pp",
    "sq",
    "tl",
]
