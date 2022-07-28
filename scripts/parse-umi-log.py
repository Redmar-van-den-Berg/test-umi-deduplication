#!/usr/bin/env python3

import argparse
import logging
import time
import os

import re
import json

def init_logger(arguments):
    # Set the log level
    # see https://docs.python.org/3.5/howto/logging.html#logging-basic-tutorial
    num_level = getattr(logging, arguments.log.upper())
    log = logging.getLogger()
    log.setLevel(num_level)

    # Console handler
    c_handler = logging.StreamHandler()
    if arguments.log.upper() == 'DEBUG':
        format_ = ' %(filename)s:%(lineno)s - %(funcName)-20s %(message)s'
    else:
        format_ = '%(levelname)-8s %(message)s'
    formatter = logging.Formatter(format_)
    c_handler.setFormatter(formatter)
    log.addHandler(c_handler)

    # File handler
    f_handler = logging.FileHandler(arguments.logfile)
    if arguments.log.upper() == 'DEBUG':
        format_ = ' %(filename)s:%(lineno)s - %(funcName)-20s %(message)s'
    else:
        format_ = '%(levelname)-8s %(message)s'
    formatter = logging.Formatter(format_)
    f_handler.setFormatter(formatter)
    log.addHandler(f_handler)

def parse(logfile):
    """ Parse umi_tools dedup log file """
    # We extract the sample name from the file name, which is fragile
    sample_name = logfile.split('.')[0]

    # The patterns we need to recognise the relevant lines
    input_reads = '.*Reads: Input Reads: \d+, Chimeric read pair: \d+$'
    output_reads = '.*Number of reads out: \d+$'
    dup_positions = '.*Total number of positions deduplicated: \d+$'
    mean_umi_per_pos = '.*Mean number of unique UMIs per position: \d+.\d+$'
    max_umi_per_pos = '.*Max. number of unique UMIs per position: \d+$'


    # Dict where we store the data
    data = dict()

    with open(logfile) as fin:
        for line in fin:
            spline = line.strip().split()

            # Get the input reads
            if re.match(input_reads, line):
                # We need to cut off the comma from the value
                data['input_reads'] = int(spline[6].replace(',',''))
                data['chimeric_read_pair'] = int(spline[10])

            # Get the output reads
            if re.match(output_reads, line):
                data['umi_tools_output'] = int(spline[7])

            # Get the duplicate positions
            if re.match(dup_positions, line):
                data['deduplicated_positions'] = int(spline[8])

            # Get the mean number of UMI's per position
            if re.match(mean_umi_per_pos, line):
                data['mean_umi_per_pos'] = float(spline[-1])

            # Get the maximum number of UMI's per position
            if re.match(max_umi_per_pos, line):
                data['max_umi_per_pos'] = int(spline[-1])

    # Calculate the percentage duplicates
    frac_unique = data['umi_tools_output'] / data['input_reads']
    frac_dup = round(1-frac_unique, 2)
    data['umi_tools_perc_duplicates'] = frac_dup * 100

    return data


def parse_stats(logfile):
    data = dict()
    with open(logfile) as fin:
        for line in fin:
            field, value = line.strip().split(': ')
            data[field] = int(value)
    # Rename some fields
    data['umi_trie_output'] = data['clusters']
    data['umi_trie_input'] = data['usable']

    frac_unique = data['clusters'] / data['usable']
    frac_dup = round(1-frac_unique, 2)
    data['umi_trie_perc_duplicates'] = frac_dup * 100
    return data

def main(args):
    data = dict()
    # Read the UMI-Tools log files
    for sample, logfile in zip(args.samples, args.umi_tools):
        results = parse(logfile)
        results['sample_name'] = sample
        data[sample] = results

    # Read the umi-trie stat files
    for sample, logfile in zip(args.samples, args.umi_trie):
        results = parse_stats(logfile)
        data[sample].update(results)


    header = ['sample_name', 'input_reads', 'umi_tools_output', 'umi_tools_perc_duplicates',
            'umi_trie_input', 'umi_trie_output', 'umi_trie_perc_duplicates']

    print(*header, sep='\t')

    for sample in data:
        print(*(data[sample][field] for field in header), sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--log', default='DEBUG', help='Level of logging info',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])

    # Some stuf we need for the logfile
    basename = os.path.basename(parser.prog)
    filename = os.path.splitext(basename)[0]
    timestamp = time.strftime('%Y%m%d_%H%M%S', time.gmtime())

    parser.add_argument('--logfile', default='/dev/null',
                        # default = '{}_{}.log'.format(filename, timestamp),
                        help='File to write the logging information to')
    parser.add_argument('--umi-tools', required=True, nargs='+',
                        help='stdout from umi_tools dedup')
    parser.add_argument('--samples', required=True, nargs='+',
                        help='Sample names, on order')
    parser.add_argument('--umi-trie', required=True, nargs='+',
                        help='stats.dat from umi-trie')

    arguments = parser.parse_args()
    assert len(arguments.umi_tools) == len(arguments.samples)
    assert len(arguments.umi_trie) == len(arguments.samples)
    init_logger(arguments)
    main(arguments)
