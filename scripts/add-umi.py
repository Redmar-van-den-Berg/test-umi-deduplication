#!/usr/bin/env python3

import argparse
import dnaio
import xopen


def main(args):
    # Open the input files
    forward = dnaio.open(args.forward, opener=xopen.xopen)
    umi = dnaio.open(args.umi, opener=xopen.xopen)
    reverse = dnaio.open(args.reverse, opener=xopen.xopen)

    # Open the output files
    fout = xopen.xopen(args.forward_out, 'wb')
    rout = xopen.xopen(args.reverse_out, 'wb')

    for f, u, r in zip(forward, umi, reverse):
        umi = u.sequence
        # If there are spaces in the header, we have to add the UMI before the
        # first space
        f_name, f_rest = f.name.split(maxsplit=1)
        r_name, r_rest = r.name.split(maxsplit=1)
        f.name = f"{f_name}_{umi} {f_rest}"
        r.name = f"{r_name}_{umi} {r_rest}"

        fout.write(f.fastq_bytes())
        rout.write(r.fastq_bytes())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('forward', help='fastq.gz file with forward reads')
    parser.add_argument('umi', help='fastq.gz file with UMI reads')
    parser.add_argument('reverse', help='fastq.gz file with reverse reads')
    parser.add_argument('--forward-out', required=True, help='Output file for forward reads')
    parser.add_argument('--reverse-out', required=True, help='Output file for reverse reads')

    arguments = parser.parse_args()
    main(arguments)
