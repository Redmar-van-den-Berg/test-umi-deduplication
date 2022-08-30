from inspect_discordant import ReadPair

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
