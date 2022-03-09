containers = {"debian": "docker://debian:latest"}
default = {"setting1": "common.smk", "setting2": "common.smk", "setting3": "common.smk"}


def get_outfile():
    return "outputfile.txt"


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
