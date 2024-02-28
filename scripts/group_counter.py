#!/usr/bin/env python3

import argparse
from typing import Iterator


def get_umis(fname: str) -> list[str]:
    umis = list()
    with open(fname) as fin:
        header = next(fin)
        for line in fin:
            umi = line.split(",")[0]
            umis.append(umi)
    return umis


def get_counts(fname: str, field: str) -> list[str]:
    counts = list()
    with open(fname) as fin:
        header = next(fin).strip("\n").split(",")
        for line in fin:
            d = {k: v for k, v in zip(header, line.strip("\n").split(","))}
            counts.append(d[field])
    return counts


def main(files: Iterator[str], include_umi: bool, field: str) -> None:
    data = dict()
    umis: list[str] = list()
    for fname in files:
        sample = fname.split("/")[0]

        if not umis:
            umis = get_umis(fname)

        data[sample] = get_counts(fname, field)

    # Set the header
    samples = list(data)
    if include_umi:
        header = ["UMI"] + samples
    else:
        header = samples
    print(*header, sep=",")

    for i in range(len(umis)):
        umi = umis[i]
        if include_umi:
            print(umi, end=",")
        for s in samples:
            print(data[s][i], end=",")
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("files", help="UMI count files", type=str, nargs="+")
    parser.add_argument(
        "--include-umi",
        default=False,
        action="store_true",
        help="Include column for UMI sequence",
    )
    parser.add_argument(
        "--field", required=True, help="Which column to include from each file"
    )

    args = parser.parse_args()

    main(args.files, args.include_umi, args.field)
