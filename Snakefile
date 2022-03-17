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
            sample=pep.sample_table.sample_name,
        ),
        bam=expand(
            "{sample}/{sample}.umi.dedup.bam", sample=pep.sample_table.sample_name
        ),
        #stats="umi-stats.tsv",


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
        umi_trie=config["umi-trie"],
    output:
        forw="{sample}/umi-trie/forward_dedup.fastq.gz",
        rev="{sample}/umi-trie/reverse_dedup.fastq.gz",
        umi="{sample}/umi-trie/umi_dedup.fastq.gz",
    log:
        "log/{sample}-umi-trie.txt",
    container:
        containers["dnaio"]
    shell:
        """
        folder=$(dirname {output.forw})
        mkdir -p $folder

        {input.umi_trie} \
            -f $folder \
            {input.forw} {input.rev} {input.umi} 2> {log}
        """


rule index_reference:
    input:
        ref=config["reference"],
    output:
        directory("gmap_index"),
    params:
        "reference",
    log:
        "log/index_reference.txt",
    container:
        containers["gsnap"]
    shell:
        """
        # Clear the existing folder
        rm -rf {output}
        mkdir {output}

        gmap_build -D {output} -d {params} {input.ref} 2>&1 > {log}
        """


rule align_vars:
    input:
        fq1=rules.concat.output.forw,
        fq2=rules.concat.output.rev,
        index=rules.index_reference.output,
    output:
        sam="{sample}/align/{sample}.sam",
    params:
        rg_sample=lambda x: "{sample}",
        db_name="reference",
    threads: 8
    log:
        "log/{sample}_align_vars.txt",
    container:
        containers["gsnap"]
    shell:
        """
        gsnap \
            --dir {input.index} \
            --db {params.db_name} \
            --batch 4 \
            --nthreads {threads} \
            --novelsplicing 1 \
            --npaths 1 \
            --quiet-if-excessive \
            --read-group-name={params.rg_sample} \
            --read-group-id={params.rg_sample} \
            --format sam \
            --gunzip {input.fq1} {input.fq2} > {output.sam} 2>{log}
        """


rule sort_bamfile:
    input:
        sam=rules.align_vars.output.sam,
    output:
        bam="{sample}/align/{sample}.bam",
        bai="{sample}/align/{sample}.bai",
    params:
        tmp=temp("tmp"),
    log:
        "log/{sample}_sort_bamfile.txt",
    container:
        containers["picard"]
    shell:
        """
        mkdir -p {params.tmp}
        picard -Xmx4G SortSam\
            I={input.sam} \
            O={output.bam} \
            SORT_ORDER=coordinate \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=true \
            TMP_DIR={params.tmp} 2>&1 > {log}
        """


rule add_umi:
    input:
        umi=rules.concat.output.umi,
        bam=rules.sort_bamfile.output.bam,
        add_umi=srcdir("scripts/bam-add-umi.py"),
    output:
        bam="{sample}/{sample}.umi.bam",
        bai="{sample}/{sample}.umi.bam.bai",
    log:
        "log/{sample}_add_umi.txt",
    container:
        containers["umi-tools"]
    shell:
        """
        python3 {input.add_umi} \
                --umi-files {input.umi} \
                --bam-file {input.bam} \
                --output-file {output.bam} 2> {log}

        python3 -c "import pysam; pysam.index('{output.bam}')"
        """


rule umi_dedup:
    input:
        bam=rules.add_umi.output.bam,
        bai=rules.add_umi.output.bai,
    output:
        bam="{sample}/{sample}.umi.dedup.bam",
    log:
        "log/{sample}_umi_dedup.log",
    container:
        containers["umi-tools"]
    shell:
        """
        mkdir -p tmp
        export TMPDIR=tmp

        umi_tools dedup \
                --stdin={input.bam} \
                --stdout={output.bam} > {log}
        """


rule parse_umi_log:
    input:
        logs=expand("log/{sample}_umi_dedup.log", sample=pep.sample_table.sample_name),
        parse_umi_log=srcdir("scripts/parse-umi-log.py"),
    output:
        "umi-stats.tsv",
    log:
        "log/parse_umi_log.txt",
    container:
        containers["umi-tools"]
    shell:
        """
        python3 {input.parse_umi_log} \
                --umi-logs {input.logs} > {output}
        """
