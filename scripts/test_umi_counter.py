import pytest
from umi_counter import umi_to_int, int_to_umi

UMIS = [
    ("ATCG", 54, 4),
    ("A", 0, 1),
    ("ATCCGATC", 13709, 8),
]


@pytest.mark.parametrize("umi, number, size", UMIS)
def test_umi_to_int(umi: str, number: int, size: int) -> None:
    assert umi_to_int(umi) == number
    assert int_to_umi(number, size) == umi
