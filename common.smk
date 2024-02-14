containers = {
    "debian": "docker://debian:latest",
    # mulled container with dnaio=0.8.1 and pysam=0.19.0
    "dnaio": "docker://quay.io/biocontainers/mulled-v2-2996a7d035117c4238b2b801e740a69df21d91e1:6b3ae5f1a97f370227e8134ba3efc0e318b288c3-0",
    "star": "docker://quay.io/biocontainers/star:2.7.10b--h9ee0642_0",
    "humid": "docker://quay.io/biocontainers/humid:1.0.2--h5f740d0_0",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.20--pyhdfd78af_0",
    "fastqc": "docker://quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1",
    "picard": "docker://quay.io/biocontainers/picard:2.27.4--hdfd78af_0",
    "umi-tools": "docker://quay.io/biocontainers/umi_tools:1.1.1--py38h0213d0e_1",
    "samtools": "docker://quay.io/biocontainers/mulled-v2-58936b48a08c4e06505165c6f560ec9460b431ea:ef260d10ee671f4c7bd8e783939839bb2e0b684e-0",
}

default = dict()
default["repeat"] = 1
default["cluster_method"] = "directional"


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("test-umi-deduplication", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Check to see if the cluster method is supported
if default["cluster_method"] not in ["directional", "maximum"]:
    raise RuntimeError(f"Unsupported cluster_method: {default['cluster_method']}")
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


def get_picard_metrics():
    metrics = list()
    for category in [
        "insert_size",
        "alignment_summary",
        "base_distribution_by_cycle",
        "quality_by_cycle",
        "quality_distribution",
    ]:
        metrics += [
            f"{sample}/align/multiple_metrics.{category}_metrics" for sample in samples
        ]

    return metrics


def get_log_files():
    """Generate paths for each of the relevant log files"""
    # HUMID log files
    humid = list()
    for category in ["stats", "clusters", "counts", "neigh"]:
        humid += [f"humid/{sample}/{category}.dat" for sample in samples]

    # UMI-Tools log files
    umi_tools = [f"log/{sample}_umi_tools.log" for sample in samples]

    # STAR log files
    STAR = [f"{sample}/align/Log.final.out" for sample in samples]

    # Picard statistics
    picard = get_picard_metrics()

    # FastQC statistics
    forward = [
        f"{sample}/concat/fastqc_{sample}/forward_fastqc.zip" for sample in samples
    ]
    reverse = [
        f"{sample}/concat/fastqc_{sample}/reverse_fastqc.zip" for sample in samples
    ]
    fastqc = forward + reverse

    return humid + umi_tools + STAR + picard + fastqc


def get_humid_after_umi_tools():
    """Generate paths for the HUMID stats file after UMI-Tools"""
    return [f"{sample}/umi-tools/humid/stats.dat" for sample in samples]


def get_umi_tools_after_humid():
    """Generate paths for the UMI-Tools log files after HUMID"""
    return [f"log/{sample}_umi_tools_after_humid.log" for sample in samples]
