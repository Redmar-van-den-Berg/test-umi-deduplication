#!/usr/bin/env python3

import argparse
from collections import Counter

NUC = {"A": "00", "C": "01", "G": "10", "T": "11"}
NUC_ = {v: k for k, v in NUC.items()}


def umi_to_int(umi: str) -> int:
    """Convert an UMI to an integer"""
    bits = list()
    for letter in umi:
        bits.append(NUC[letter])
    return int("".join(bits), 2)


def int_to_umi(number: int, size: int) -> str:
    """Convert number to an UMI of the specified size"""
    bits = format(number, "b").zfill(size * 2)
    umi = list()
    for i in range(0, len(bits), 2):
        bit = bits[i : i + 2]
        umi.append(NUC_[bit])

    return "".join(umi)


def extract_umi(read):
    """Extract the UMI from the read header"""
    return read.name.split(" ")[0].split(":")[-1]


def fetch_umis(fname):
    with dnaio.open(fname, opener=xopen.xopen, mode="r") as fin:
        for read in fin:
            yield extract_umi(read)


def fetch_umis_as_int(fname):
    for umi in fetch_umis(fname):
        try:
            yield umi_to_int(umi)
        except KeyError:
            continue


def write_alphabetic(counts, umi_size, fname):
    total = counts.total()
    cumulative = 0
    with open(fname, "wt") as fout:
        print("UMI", "count", "fraction", "cumulative", "cumulative_frac", file=fout, sep=",")
        for i in range(4**umi_size):
            umi = int_to_umi(i, umi_size)
            count = counts[i]
            fraction = count / total

            # Cumulative
            cumulative += count
            cumulative_frac = cumulative / total

            print(umi, count, fraction, cumulative, cumulative_frac, file=fout, sep=",")


def write_descending(counts, umi_size, fname):
    total = counts.total()
    cumulative = 0
    with open(fname, "wt") as fout:
        print("UMI", "count", "fraction", "cumulative", "cumulative_frac", file=fout, sep=",")
        for i, count in counts.most_common():
            umi = int_to_umi(i, umi_size)
            count = counts[i]
            fraction = count / total

            # Cumulative
            cumulative += count
            cumulative_frac = cumulative / total

            print(umi, count, fraction, cumulative, cumulative_frac, file=fout, sep=",")


def calculate_bias(counts, umi_size):
    """Calculate the bias from a perfectly even distribution of UMI's"""
    total = counts.total()
    possible_umi = 4**umi_size
    target_count = total / possible_umi
    bias = 0

    for _, count in counts.items():
        bias += abs(count - target_count)

    return bias, bias/total


def main(fname, umi_size, alphabetic_file, descending_file):
    # Initialise the counter with all possible UMI's set to zero
    counts = Counter({x:0 for x in range(4**umi_size)})
    counts.update(fetch_umis_as_int(fname))

    write_alphabetic(counts, umi_size, alphabetic_file)
    write_descending(counts, umi_size, descending_file)

    bias_umi, bias_umi_fraction = calculate_bias(counts, umi_size)

    print(f"{bias_umi=}")
    print(f"{bias_umi_fraction=}")


if __name__ == "__main__":
    import dnaio
    import xopen
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Count UMIs in FastQ header",
    )

    parser.add_argument("fastq", help="Fastq input file")
    parser.add_argument("--umi-size", type=int, default=8, help="Size of the UMI")
    parser.add_argument(
        "--alphabetic",
        required=True,
        type=str,
        help="Write UMI counts to this file in alphabetic order",
    )
    parser.add_argument(
        "--descending",
        required=True,
        type=str,
        help="Write UMI counts to this file in descending order",
    )
    args = parser.parse_args()
    main(args.fastq, args.umi_size, args.alphabetic, args.descending)
