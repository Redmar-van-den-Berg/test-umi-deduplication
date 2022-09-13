#!/usr/bin/env python3

from dataclasses import dataclass, field
from collections import defaultdict

import argparse
import functools
import pysam
import sys
import json

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
    for cluster in find_matches(clusters, item, compare):
        new_cluster += cluster
        clusters.remove(cluster)

    # Add item itself to the new cluster
    new_cluster.append(item)
    
    # Add the new cluster to the list of clusters
    clusters.append(new_cluster)


def explain_discordance(cluster):
    """ Attempt to explain the discordance between the reads in the cluster """
    if len(cluster) == 2:
        pair1 = cluster[0]
        pair2 = cluster[1]

        # If both tools simply picked a different read to mark as duplicate
        if (pair1.trie_filt and pair2.tool_filt) or (pair1.tool_filt and pair2.trie_filt):
            return 'Alternative read'

    if len(cluster) == 1:
        pair = cluster[0]
        if pair.trie_filt:
            return 'Single read, marked by umi-trie'
        if pair.tool_filt:
            return 'Single read, marked by umi-tools'


def write_unexplained(infile, outfile, unexplained):
    samfile = pysam.AlignmentFile(infile)
    outfile = pysam.AlignmentFile(outfile, "wb", template=samfile)

    for read in samfile.fetch():
        if read.query_name in unexplained:
            outfile.write(read)

    samfile.close()
    outfile.close()
    pysam.index(outfile)

def main(args):
    samfile = pysam.AlignmentFile(args.bam)

    readpairs = defaultdict(list)

    # Result of explanations for discordant reads
    results = defaultdict(int)
    # Keep track of all reads that we could not explain
    unexplained = list()

    for i, read in enumerate(samfile.fetch(),1):
        if i % 1000 == 0:
            print(f'Parsed {i} reads from {args.bam}', file=sys.stderr)
        if i % 20000 == 0:
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
        clusters = list()
        comp = functools.partial(match_reads, distance=1)
        # Merge reads together based on hamming distance
        for read in reads:
            add_and_merge(clusters, read, comp)

        for cluster in clusters:
            if (exp := explain_discordance(cluster)):
                results[exp] +=2
            else:
                results['unexplained'] += len(cluster)
                unexplained += cluster
                print(pos)
                print(*cluster, sep='\n')
                print()


    # Put the readnames for unxplained reads in a set
    unex = {read.name for read in unexplained}

    # Write unexplained reads to bam output
    print(f'Writing {len(unex)} reads to {args.bamout}', file=sys.stderr)
    write_unexplained(args.bam, args.bamout, unex)

    print(json.dumps(results,indent=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('bam', help='Input bamfile')
    parser.add_argument('bamout', help='Output bamfile')

    args = parser.parse_args()
    main(args)
