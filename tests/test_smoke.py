import pytest

import scanpex as scx


def test_smoke():
    assert scx.__version__ is not None

    from scanpex import ft, ml, pl, pp, sq, tl
