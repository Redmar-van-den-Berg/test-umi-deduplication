include: "common.smk"


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("test-umi-deduplication", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Set the updated dict as the configuration for the pipeline
config = default


rule all:
    input:
        concat=expand(
            "{sample}/umi-trie/forward_dedup.fastq.gz",
            sample=pep.sample_table["sample_name"],
        ),


rule concat:
    """Concatentate the input fastq files"""
    input:
        forw=get_forward,
        rev=get_reverse,
        umi=get_umi,
    output:
        forw="{sample}/concat/forward.fastq.gz",
        rev="{sample}/concat/reverse.fastq.gz",
        umi="{sample}/concat/umi.fastq.gz",
    log:
        "log/{sample}_concat.txt",
    container:
        containers["debian"]
    shell:
        """
        mkdir -p $(dirname {output.forw})

        cp {input.forw} {output.forw} || cat {input.forw} > {output.forw}
        cp {input.rev} {output.rev} || cat {input.rev} > {output.rev}
        cp {input.umi} {output.umi} || cat {input.umi} > {output.umi}
        """


rule umi_trie:
    """Run umi-trie on the fastq files"""
    input:
        forw=rules.concat.output.forw,
        rev=rules.concat.output.rev,
        umi=rules.concat.output.umi,
        umi_trie=srcdir("bin/umi-trie"),
    output:
        forw="{sample}/umi-trie/forward_dedup.fastq.gz",
        rev="{sample}/umi-trie/reverse_dedup.fastq.gz",
        umi="{sample}/umi-trie/umi_dedup.fastq.gz",
    log:
        "log/{sample}-umi-trie.txt",
    container:
        containers["debian"]
    shell:
        """
        folder=$(dirname {output.forw})
        mkdir -p $folder

        {input.umi_trie} \
            -f $folder \
            {input.forw} {input.rev} {input.umi} 2> {log}
        """
