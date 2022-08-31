import pytest

from inspect_discordant import add_and_merge, find_matches, ReadPair, match_word,  match_reads, explain_discordance

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

def test_match_read():
    assert match_reads(readpair1, readpair2)

def test_add_and_merge_empty():
    clusters = list()
    item = 'AA'
    add_and_merge(clusters, item, match_word)
    assert clusters == [['AA']]

def test_add_and_merge_match():
    clusters = [['AB']]
    item = 'AA'
    add_and_merge(clusters, item, match_word)
    assert clusters == [['AB', 'AA']]

def test_add_and_merge_mismatch():
    clusters = [['BB']]
    item = 'AA'
    add_and_merge(clusters, item, match_word)
    assert clusters == [['BB'], ['AA']]

def test_add_and_merge_bridge():
    """ Test if two clusters are merged if item bridges them """
    clusters = [['AA'], ['BB']]
    item = 'AB'
    add_and_merge(clusters, item, match_word)
    assert clusters == [['AA', 'BB', 'AB']]

def test_find_matches_none():
    clusters = list()
    item = 'AA'
    find_matches(clusters, item, match_word) == []

def test_find_matches_one():
    clusters = [['AA']]
    item = 'AB'
    assert find_matches(clusters, item, match_word) == [['AA']]

def test_find_matches_two():
    clusters = [['AA'], ['BB']]
    item = 'AB'
    assert find_matches(clusters, item, match_word) == [['AA'], ['BB']]

def test_explain_discordance():
    cluster = [readpair1, readpair2]
    assert explain_discordance(cluster) == 'Alternative read'
