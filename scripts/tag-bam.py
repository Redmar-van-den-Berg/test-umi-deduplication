#!/usr/bin/env python3

import argparse
import logging
import time
import os
import sys

import dnaio
import pysam
import xopen

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


def get_bam_reads(filename):
    """Extract the read names from a bam file"""
    reads = set()
    samfile = pysam.AlignmentFile(filename, 'rb')

    for read in samfile:
        header = read.query_name.split('_')[0]

        # Print example header the first time
        if not reads:
            logging.info(f'Storing first bam header: "{header}"')

        reads.add(header)

    return reads


def get_fastq_reads(filename):
    reads = set()

    fin = dnaio.open(filename, opener=xopen.xopen, mode='r')

    for read in fin:
        header = read.name.split(' ')[0]

        if not reads:
            logging.info(f'Storing the first fastq header: "{header}"')
        reads.add(header)

    return reads

def main(args):
    bam_reads = get_bam_reads(args.bam_reads)
    fastq_reads = get_fastq_reads(args.fastq_reads)

    # Open the bamfile that we want to tag
    bam_in = pysam.AlignmentFile(args.bam_in, 'rb')
    bam_out = pysam.AlignmentFile(args.bam_out, 'wb', template=bam_in)


    # Colors for the tagged reads
    colors = {
            0: '0,136,86', # Unfiltered, green
            1: '0,103,165', # UMI-trie filtered, blue
            2: '135,86,146', # UMI-Tools filtered, purple
            3: '190,0,50' # Filtered by both, red 
    }
    first = True
    for i,read in enumerate(bam_in.fetch(), start=1):
        readname = read.query_name.split('_')[0]

        if first:
            logging.info(f'First header from bam file: "{readname}"')
            first=False

        filtered = 0
        # If filtered by UMI-Tools
        if readname not in bam_reads:
            filtered += 1
        # If filtered by umi-trie
        if readname not in fastq_reads:
            filtered += 2

        c = colors[filtered]
        read.tags +=[('XU:i', filtered), ('YC:Z', c)]
        bam_out.write(read)

        if i%100000 == 0:
            print(f'Written {i} reads to {args.bam_out}\r')
    logging.info(f'Written {i} reads to {args.bam_out}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Set tag based on wether the tool filtered a read. UMITOOLS 0 means umitools did not filter this read'
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

    parser.add_argument('--bam-in', required=True, help='BAM file to tag')
    parser.add_argument('--fastq-reads', required=True, help='FASTQ file to get the reads from')
    parser.add_argument('--bam-reads', required=True, help='BAM file to get the reads from')
    parser.add_argument('--bam-out', required=True, help='Tagged bam file')

    arguments = parser.parse_args()
    init_logger(arguments)
    main(arguments)
