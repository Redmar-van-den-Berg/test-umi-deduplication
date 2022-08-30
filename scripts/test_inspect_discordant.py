import pytest

from inspect_discordant import ReadPair, match_word

readpair1 = ReadPair('name_ATCG', (1,2), 1, 'ATCG')
readpair2 = ReadPair('name2_ATCG', (1,2), 2, 'ATCG')
def test_true():
    assert True

def test_trie_filt():
    """ Test field for umi-trie filter """
    assert readpair1.trie_filt

def test_tool_filt():
    """ Test field for umi-tools filter """
    assert readpair2.tool_filt

def test_match_word():
    """ Test if two words match. They must be the same length """
    assert match_word('', '')
    assert match_word('', '', 0)
    assert match_word('A', 'A')
    assert match_word('A', 'B')
    assert not match_word('A', 'B', 0)
    with pytest.raises(RuntimeError):
        match_word('AA', 'B')
