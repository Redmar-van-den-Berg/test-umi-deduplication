containers = {
    "debian": "docker://debian:latest",
    "dnaio": "docker://quay.io/biocontainers/dnaio:0.7.1--py39hbf8eff0_1",
    "gsnap": "docker://quay.io/biocontainers/gmap:2014.12.23--pl5.22.0_4",
    "picard": "docker://quay.io/biocontainers/picard:2.20.5--0",
    "umi-tools": "docker://quay.io/biocontainers/umi_tools:1.1.1--py38h0213d0e_1",
}

default = {"umi-trie": srcdir("bin/umi-trie")}


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("test-umi-deduplication", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Set the updated dict as the configuration for the pipeline
config = default

# Make sample names easily accessible
samples = list(pep.sample_table.sample_name)


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
