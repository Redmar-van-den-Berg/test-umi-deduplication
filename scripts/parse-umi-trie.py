#!/usr/bin/env python3

import argparse
import re
import json

def parse(fin):
    data = dict()
    for line in fin:
        field, value = line.strip().split(': ')
        data[field] = int(value)
    # Rename some fields
    data['umi_trie_output'] = data['clusters']
    data['umi_trie_input'] = data['usable']

    frac_unique = data['clusters'] / data['usable']
    frac_dup = round((1-frac_unique)*100, 2)
    data['Percentage of duplicate reads'] = frac_dup
    return data

def main(args):
    # Read the umi-trie stat files
    with open(args.umitrie) as fin:
        print(json.dumps(parse(fin), indent=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('umitrie', help='stats.dat from umi-trie')

    arguments = parser.parse_args()
    main(arguments)
