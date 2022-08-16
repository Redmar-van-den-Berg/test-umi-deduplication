#!/usr/bin/env python3

import argparse
import re
import json

def split_value(value):
    field, value = value.split(': ')
    try:
        value = int(value)
    except ValueError:
        try:
            value = float(value)
        except ValueError:
            pass
    return field, value

def parse_line(line):
    line = line.strip().split('INFO ')[1]
    return split_value(line)

def parse(fin):
    """ Parse umi_tools dedup log file """
    # The patterns we need to recognise the relevant lines
    input_reads = '.*Reads: Input Reads: \d+, Chimeric read pair: \d+$'
    output_reads = '.*Number of reads out: \d+$'
    dup_positions = '.*Total number of positions deduplicated: \d+$'
    mean_umi_per_pos = '.*Mean number of unique UMIs per position: \d+.\d+$'
    max_umi_per_pos = '.*Max. number of unique UMIs per position: \d+$'


    # Dict where we store the data
    data = dict()

    for line in fin:
        # Get the input reads
        if re.match(input_reads, line):
            # We need to cut off the comma from the value
            # Cut off the INFO preabmle
            line = line.split('INFO Reads: ')[1]

            # There are two key/value pairs in this line
            input_, chimeric = line.split(', ')

            # First pair, input read
            field, value = split_value(input_)
            data[field] = value

            # Second pair, chimeric reads
            field, value = split_value(chimeric)
            data[field] = value

        # Get the output reads
        if re.match(output_reads, line):
            field, value = parse_line(line)
            data[field] = value

        # Get the duplicate positions
        if re.match(dup_positions, line):
            field, value = parse_line(line)
            data[field] = value

        # Get the mean number of UMI's per position
        if re.match(mean_umi_per_pos, line):
            field, value = parse_line(line)
            data[field] = value

        # Get the maximum number of UMI's per position
        if re.match(max_umi_per_pos, line):
            field, value = parse_line(line)
            data[field] = value

    # Calculate the percentage duplicates
    frac_unique = data['Number of reads out'] / data['Input Reads']
    frac_dup = round((1-frac_unique)*100, 2)
    data['Percentage of duplicate reads'] = frac_dup

    return data


def main(args):
    data = dict()
    # Read the UMI-Tools log files
    with open(args.umitools) as fin:
        print(json.dumps(parse(fin), indent=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('umitools', help='stdout from umi_tools dedup')

    arguments = parser.parse_args()
    main(arguments)
