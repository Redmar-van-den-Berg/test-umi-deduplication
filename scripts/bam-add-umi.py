#!/usr/bin/env python3

import argparse
import logging
import time
import os

import gzip
import pysam
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

def extract_read_id(header):
    """ Extract the read ID from the header

    The header is of format:
    @header optional more stuf
    """
    return header[1:].split(' ')[0]

def read_umis(filename):
    """ Read umi's from filename and return as a dict """
    umis = dict()
    with gzip.open(filename, 'rt') as fin:
        while True:
            try:
                read_id = extract_read_id(next(fin))
            except StopIteration:
                break
            umi = next(fin).strip()
            umis[read_id] = umi

            # Skip second header and qual
            next(fin)
            next(fin)

            # Print status message
            if len(umis) % 1000000 == 0:
                logging.info(f"Parsed {len(umis)} UMI's from {filename}")

    return umis

def write_with_umis(umis, bamfile, outfile):
    """ Add the umis to the bamfiles and write to outfile """
    print(f"Working with {len(umis)} umis")
    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        with pysam.AlignmentFile(outfile, 'wb', template=samfile) as output:
            for i, record in enumerate(samfile):
                # Add the umi to the read name
                try:
                    umi = umis[record.qname]
                except KeyError:
                    logging.error(f'Unkown header "{record.qname}", skipping')
                    continue
                record.qname += f'_{umi}'
                output.write(record)

                # Print status message
                if i % 1000000 == 0:
                    logging.info(f"Wrote {i} reads to {outfile}")

def main(args):

    # Read the UMIs from the specified UMI files if we didn't get a json
    if not args.umi_json:
        umis = dict()
        for filename in args.umi_files:
            umis.update(read_umis(filename))
        #with open('umi.json', 'wt') as fout:
        #    print(json.dumps(umis, indent=2), file=fout)
    else:
        with open(args.umi_json) as fin:
            umis = json.load(fin)

    write_with_umis(umis, args.bam_file, args.output_file)

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
    parser.add_argument('--umi-files', required=True, nargs='+',
                        help='fastq.gz file(s) with UMIs, as provided by '
                        'GenomeScan')
    parser.add_argument('--bam-file', required=True,
                        help='Aligned and sorted bam file')
    parser.add_argument('--output-file', required=True,
                        help='Bam file to write the output to')
    parser.add_argument('--umi-json', required=False,
                        help='JSON umi file')

    arguments = parser.parse_args()
    init_logger(arguments)
    main(arguments)
