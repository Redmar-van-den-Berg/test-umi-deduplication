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
    return read.seq[:8] + pair.seq[:8] + umi


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

