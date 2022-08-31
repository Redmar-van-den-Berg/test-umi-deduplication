#!/usr/bin/env python3

from dataclasses import dataclass, field
from collections import defaultdict

import pysam
import sys

@dataclass
class ReadPair:
    name: str
    pos: tuple
    XU: int
    word: str = field(compare=False)
    trie_filt: bool = field(init=False)
    tool_filt: bool = field(init=False)

    def __post_init__(self):
        self.trie_filt = self.XU % 2 == 1
        self.tool_filt = self.XU >= 2

def get_chr(read, samfile):
    return samfile.get_reference_name(read.reference_id)

def make_pos(read, pair, samfile):
    """ Return a tuple of the positions of read and pair

    """
    assert get_chr(read, samfile) == get_chr(pair, samfile)
    return (read.pos, pair.pos)


def make_word(read, pair):
    """ Create a word, similar to umi-trie """
    umi = read.query_name.split('_')[1]
    return umi + read.seq[:8] + pair.seq[:8]


def make_readpair(read, samfile):
    """ Make a readpair out of a read

    Assumes that the mate is properly paired
    """
    assert read.is_proper_pair

    pair = samfile.mate(read)
    pos = make_pos(read, pair, samfile)
    word = make_word(read, pair)

    return ReadPair(read.query_name, pos, read.get_tag('XU'), word)

def mate_unmapped(read, samfile):
    """ Is the mate of read one of the various types to be unmapped"""

    # The readpair is not proper, e.g. mapped to different chromosomes
    if not read.is_proper_pair:
        return True

    # The mate cannot be found (testing only)
    try:
        samfile.mate(read)
    except ValueError:
        return True

    return False


def is_discordant(read):
    """ Determine if a read is discordant between UMI-Tools and umi-trie """
    xu = read.get_tag('XU')
    return xu == 1 or xu == 2


def match_word(word1, word2, distance=1):
    """ Does word1 match word2, given the specified distance """
    # Words must be the same size
    if len(word1) != len(word2):
        raise RuntimeError(f'Unequal word length: "{word1}" and "{word2}"')

    d = 0
    for i, j in zip(word1, word2):
        if i != j:
            d+=1
        if d > distance:
            return False
    return True


def match_reads(read1, read2, distance=1):
    return match_word(read1.word, read2.word, distance=distance)


def find_matches(clusters, item, compare):
    """ Determine which clusters match item, using compare """
    matches = list()

    for cluster in clusters:
        if any(compare(item, i) for i in cluster):
            matches.append(cluster)
            continue
    return matches


def add_and_merge(clusters, item, compare):
    """ Add item to one of the clusters, or add it to a new cluster

    Use compare to determine if item fits with any of the clusters.

    If the addition of item bridges multiple clusters, merge them

    This function modifies clusters in place
    """
    # We always create a new cluster. It will either contain only item, or any
    # existing clusters that match item, as well as item itself
    new_cluster = list()

    # Find and merge all existing cluster that match item
    matches = find_matches(clusters, item, compare)
    for cluster in matches:
        new_cluster += cluster
        clusters.remove(cluster)

    # Add item itself to the new cluster
    new_cluster.append(item)
    
    # Add the new cluster to the list of clusters
    clusters.append(new_cluster)

if __name__ == '__main__':
    fname = sys.argv[1]

    samfile = pysam.AlignmentFile(fname)

    readpairs = defaultdict(list)

    # Result of explanations for discordant reads
    results = defaultdict(int)

    for i, read in enumerate(samfile.fetch(),1):
        if i % 1000 == 0:
            print(f'Parsed {i} reads from {fname}', file=sys.stderr)
            break

        # We only act on read1 (since we fetch the mate)
        if not read.is_read1:
            continue

        # Skip reads that are not discordant between the two tools
        if not is_discordant(read):
            continue

        if mate_unmapped(read, samfile):
            results['mate unmapped']+=1
            continue

        # Convert into a readpair object
        readpair = make_readpair(read, samfile)

        # Store for later parsing
        readpairs[readpair.pos].append(readpair)

    for pos, reads in readpairs.items():
        print(pos)
        print(*reads, sep='\n')
        print()

