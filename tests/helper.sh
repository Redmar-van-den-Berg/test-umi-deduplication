#!/usr/bin/env bash

set -euxo pipefail

function setup {
  unxz -k tests/data/star/Genome.xz
  unxz -k tests/data/reference/hamlet-ref.fa.xz
}

function cleanup {
  rm -f tests/data/star/Genome
  rm -f tests/data/reference/hamlet-ref.fa
}

trap cleanup EXIT

setup
