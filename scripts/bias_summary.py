#!/usr/bin/env python3

import argparse

def read_bias(fname):
    data = dict()
    with open(fname) as fin:
        for line in fin:
            field, value = line.strip('\n').split('=')
            data[field] = value
    return data

def main(files, names):
    assert len(files) == len(names)

    header = ["bias_umi", "bias_umi_fraction"]
    print(*header, sep=',')

    for fname, name in zip(files, names):
        d = read_bias(fname)
        print(*(d[field] for field in header), sep=',')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Count UMIs in FastQ header",
    )

    parser.add_argument("--bias", nargs='+', help="Files that contain the bias")
    parser.add_argument("--samples", nargs='+', help="Sample names, in the same order")

    args = parser.parse_args()
    main(args.bias, args.samples)
