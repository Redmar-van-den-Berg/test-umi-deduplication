#!/usr/bin/env python3

import argparse
import logging
import time
import os

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
        reads.add(read.query_name.split('_')[0])

    return reads


def main(args):
    reads = get_bam_reads(args.bam)
    log_n_reads=1000000

    for infile, outfile in zip(args.fastq_in, args.fastq_out):
        fin = dnaio.open(infile, opener=xopen.xopen, mode='r')
        fout = xopen.xopen(outfile, mode='wb')

        counter = 0
        # Keep track of the reads we have seen for this file
        file_reads = reads.copy()

        for read in fin:
            # Print status message
            counter+=1
            if counter%log_n_reads== 0:
                logging.info(f'Parsed {counter/log_n_reads}M reads from {infile}\r')

            # Take the first part of the header and cut of the @
            header_name = read.name.split(' ')[0]

            if header_name in file_reads:
                fout.write(read.fastq_bytes())
                # Remove the header once we have seen it, to determine if we missed any once we are done
                file_reads.remove(header_name)


        if file_reads:
            logging.error(f'{len(reads)} reads could not be found\n')
            logging.error(f'Example: {reads.pop()}')
            exit(1)
        else:
            logging.info(f'Parsed {counter/log_n_reads}M reads from {infile}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Intersect fastq files with a bam file'
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

    parser.add_argument('--fastq-in', required=True, nargs='+', help='Fastq input files')
    parser.add_argument('--fastq-out', required=True, nargs='+', help='fastq output files')
    parser.add_argument('--bam', required=True, help='BAM input file')
    arguments = parser.parse_args()
    assert len(arguments.fastq_in) == len(arguments.fastq_out)
    main(arguments)
