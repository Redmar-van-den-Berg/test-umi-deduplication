containers = {
    "debian": "docker://debian:latest",
    "dnaio": "docker://quay.io/biocontainers/dnaio:0.7.1--py39hbf8eff0_1",
    "gsnap": "docker://quay.io/biocontainers/gmap:2014.12.23--pl5.22.0_4",
    "picard": "docker://quay.io/biocontainers/picard:2.20.5--0",
}

default = {"umi-trie": srcdir("bin/umi-trie")}


def get_fastq(wildcards, column):
    fastq = pep.sample_table.loc[wildcards.sample, column]

    # If a single fastq file is specified, forward will be a string
    if isinstance(fastq, str):
        return [fastq]
    # If multiple fastq files were specified, forward will be a list
    else:
        return fastq


def get_forward(wildcards):
    return get_fastq(wildcards, "forward")


def get_reverse(wildcards):
    return get_fastq(wildcards, "reverse")


def get_umi(wildcards):
    return get_fastq(wildcards, "umi")
