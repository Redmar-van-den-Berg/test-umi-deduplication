include: "common.smk"


rule all:
    input:
        concat=expand(
            "humid/{sample}/forward_dedup.fastq.gz",
            sample=samples,
        ),
        bam=expand("{sample}/{sample}.umi.dedup.bam", sample=samples),
        stats="umi-stats.tsv",
        fastq=expand("{sample}/umi-tools/forward.fastq.gz", sample=samples),
        trie_after_umi_tools=expand(
            "{sample}/umi-tools/umi-trie/forward_dedup.fastq.gz", sample=samples
        ),
        umi_tools_after_trie=expand(
            "{sample}/umi-trie/umi-tools/{sample}.umi.dedup.bam", sample=samples
        ),
        multiqc="multiqc_report.html",
        benchmarks="benchmarks/s.tsv",


rule concat:
    """Concatentate the input fastq files"""
    input:
        forw=get_forward,
        rev=get_reverse,
        umi=get_umi,
    output:
        forw=temp("{sample}/concat/forward.fastq.gz"),
        rev=temp("{sample}/concat/reverse.fastq.gz"),
        umi=temp("{sample}/concat/umi.fastq.gz"),
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


rule humid:
    """Run umi-trie on the fastq files"""
    input:
        forw=rules.concat.output.forw,
        rev=rules.concat.output.rev,
        umi=rules.concat.output.umi,
    params:
        cluster_method="-x" if config["cluster_method"] == "maximum" else "",
    output:
        forw="humid/{sample}/forward_dedup.fastq.gz",
        rev="humid/{sample}/reverse_dedup.fastq.gz",
        umi="humid/{sample}/umi_dedup.fastq.gz",
        stats="humid/{sample}/stats.dat",
    log:
        "log/{sample}-humid.txt",
    benchmark:
        repeat("benchmarks/humid_{sample}.tsv", config["repeat"])
    container:
        containers["humid"]
    shell:
        """
        folder=$(dirname {output.forw})
        mkdir -p $folder

        humid \
            {params.cluster_method} \
            -d $folder \
            -s \
            {input.forw} {input.rev} {input.umi} 2> {log}
        """


rule add_umi:
    input:
        forw=rules.concat.output.forw,
        umi=rules.concat.output.umi,
        rev=rules.concat.output.rev,
        add_umi=srcdir("scripts/add-umi.py"),
    output:
        forw="{sample}/{sample}.umi.forward.fastq.gz",
        rev="{sample}/{sample}.umi.reverse.fastq.gz",
    log:
        "log/{sample}_add_umi.txt",
    container:
        containers["dnaio"]
    shell:
        """
        python3 {input.add_umi} \
                {input.forw} \
                {input.umi} \
                {input.rev} \
                --forward-out {output.forw} \
                --reverse-out {output.rev} 2> {log}
        """


rule align_vars:
    input:
        fq1=rules.add_umi.output.forw,
        fq2=rules.add_umi.output.rev,
        index=config["star_index"],
        gtf=config["gtf"],
    output:
        bam="{sample}/align/Aligned.sortedByCoord.out.bam",
    params:
        rg_sample=lambda wildcards: wildcards.sample,
        chim_segment=20,
        min_intron_size=50,
        alignInsertionFlush="Right",
        twopassMode="Basic",
    log:
        main="{sample}/snv-indels/Log.out",
        progress="{sample}/snv-indels/Log.progress.out",
        final="{sample}/snv-indels/Log.final.out",
    benchmark:
        repeat("benchmarks/STAR_{sample}.tsv", config["repeat"])
    threads: 8
    container:
        containers["star"]
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --sjdbGTFfile {input.gtf} \
            --readFilesCommand zcat \
            --outSAMattrRGline "ID:{params.rg_sample}" "SM:{params.rg_sample}" \
            --outFileNamePrefix $(dirname {output.bam})/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --alignIntronMin {params.min_intron_size} \
            --alignInsertionFlush {params.alignInsertionFlush} \
            --twopassMode {params.twopassMode} \
            --chimOutType WithinBAM \
            --chimSegmentMin {params.chim_segment} \
            --quantMode GeneCounts \
            --readFilesIn {input.fq1:q} {input.fq2:q}
        """


rule index_bamfile:
    input:
        bam=rules.align_vars.output.bam,
    output:
        bai="{sample}/align/Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/index_bamfile.{sample}.txt",
    container:
        containers["samtools"]
    shell:
        """
        samtools index {input.bam} 2> {log}
        """


rule umi_tools:
    input:
        bam=rules.align_vars.output.bam,
        bai=rules.index_bamfile.output.bai,
    params:
        cluster_method="--method cluster"
        if config["cluster_method"] == "maximum"
        else "--method directional",
    output:
        bam="{sample}/{sample}.umi.dedup.bam",
    log:
        "log/{sample}_umi_dedup.log",
    benchmark:
        repeat("benchmarks/umi_tools_{sample}.tsv", config["repeat"])
    container:
        containers["umi-tools"]
    shell:
        """
        mkdir -p tmp
        export TMPDIR=tmp

        umi_tools dedup \
                {params.cluster_method} \
                --stdin={input.bam} \
                --stdout={output.bam} > {log}
        """


rule parse_umi_log:
    input:
        umi_tools=expand("log/{sample}_umi_dedup.log", sample=samples),
        umi_trie=expand("humid/{sample}/stats.dat", sample=samples),
        parse_umi_log=srcdir("scripts/parse-umi-log.py"),
    params:
        samples=samples,
    output:
        "umi-stats.tsv",
    log:
        "log/parse_umi_log.txt",
    container:
        containers["umi-tools"]
    shell:
        """
        python3 {input.parse_umi_log} \
                --samples {params.samples} \
                --umi-trie {input.umi_trie} \
                --umi-tools {input.umi_tools} > {output} 2> {log}
        """


rule bam_to_fastq:
    input:
        forw=rules.concat.output.forw,
        rev=rules.concat.output.rev,
        umi=rules.concat.output.umi,
        bam=rules.umi_tools.output.bam,
        src=srcdir("scripts/fastq-from-bam.py"),
    output:
        forw="{sample}/umi-tools/forward.fastq.gz",
        rev="{sample}/umi-tools/reverse.fastq.gz",
        umi="{sample}/umi-tools/umi.fastq.gz",
    log:
        "log/{sample}.bam_to_fastq.txt",
    container:
        containers["dnaio"]
    shell:
        """
        python3 {input.src} \
            --fastq-in {input.forw} {input.rev} {input.umi} \
            --fastq-out {output.forw} {output.rev} {output.umi} \
            --bam {input.bam} 2> {log}
        """


# Run umi-trie on the output of umi-tools
use rule humid as humid_after_umi_tools with:
    input:
        forw=rules.bam_to_fastq.output.forw,
        rev=rules.bam_to_fastq.output.rev,
        umi=rules.bam_to_fastq.output.umi,
    output:
        forw="{sample}/umi-tools/umi-trie/forward_dedup.fastq.gz",
        rev="{sample}/umi-tools/umi-trie/reverse_dedup.fastq.gz",
        umi="{sample}/umi-tools/umi-trie/umi_dedup.fastq.gz",
        stats="{sample}/umi-tools/umi-trie/stats.dat",
    log:
        "log/{sample}.humid_after_umi_tools.txt",
    benchmark:
        repeat("benchmarks/humid_after_umi_tools_{sample}.tsv", config["repeat"])


# Add the UMI's to the read name after running UMI-trie
use rule add_umi as add_umi_after_umi_trie with:
    input:
        forw=rules.humid.output.forw,
        rev=rules.humid.output.rev,
        umi=rules.humid.output.umi,
        add_umi=srcdir("scripts/add-umi.py"),
    output:
        forw="{sample}/umi-trie/umi-tools/{sample}.umi.forward.fastq.gz",
        rev="{sample}/umi-trie/umi-tools/{sample}.umi.reverse.fastq.gz",
    log:
        "log/{sample}.add_umi_after_umi_trie.txt",


# Align the reads to the reference
use rule align_vars as align_after_umi_trie with:
    input:
        fq1=rules.add_umi_after_umi_trie.output.forw,
        fq2=rules.add_umi_after_umi_trie.output.rev,
        index=config["star_index"],
        gtf=config["gtf"],
    output:
        bam="{sample}/umi-trie/umi-tools/align/Aligned.sortedByCoord.out.bam",
    log:
        "log/{sample}.align_after_umi_trie.txt",
    benchmark:
        repeat("benchmarks/STAR_after_umi_trie_{sample}.tsv", config["repeat"])


use rule index_bamfile as index_after_umi_trie with:
    input:
        bam=rules.align_after_umi_trie.output.bam,
    output:
        bai="{sample}/umi-trie/umi-tools/align/Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/{sample}.index_after_umi_trie.txt",


# Run umi-tools
use rule umi_tools as umi_tools_after_umi_trie with:
    input:
        bam=rules.align_after_umi_trie.output.bam,
        bai=rules.index_after_umi_trie.output.bai,
    output:
        bam="{sample}/umi-trie/umi-tools/{sample}.umi.dedup.bam",
    log:
        "log/{sample}_umi_dedup_after_umi_trie.log",
    benchmark:
        repeat("benchmarks/umi_tools_after_umi_trie_{sample}.tsv", config["repeat"])


# Run MultiQC on HUMID output
rule multiqc:
    input:
        stats=get_log_files(),
        config=srcdir("cfg/multiqc.yml"),
    output:
        html="multiqc_report.html",
    log:
        "log/multiqc.txt",
    container:
        containers["multiqc"]
    shell:
        """
        rm -f multiqc_filelist.txt

        for fname in {input.stats}; do
            echo $fname >> multiqc_filelist.txt
        done

        multiqc \
        --force \
        --file-list multiqc_filelist.txt \
        --config {input.config} 2> {log}
        """


rule gather_benchmarks:
    input:
        humid=expand("benchmarks/humid_{sample}.tsv", sample=samples),
        star=expand("benchmarks/STAR_{sample}.tsv", sample=samples),
        umi_tools=expand("benchmarks/umi_tools_{sample}.tsv", sample=samples),
        script=srcdir("scripts/parse-benchmark.py"),
    params:
        samples=samples,
    output:
        seconds="benchmarks/s.tsv",
        cpu_time="benchmarks/cpu_time.tsv",
        max_rss="benchmarks/max_rss.tsv",
    log:
        "log/gather_benchmarks.txt",
    container:
        containers["dnaio"]
    shell:
        """
        for column in s cpu_time max_rss; do
          python3 {input.script} \
              --samples {params.samples} \
              --tools humid STAR umi_tools \
              --column ${{column}} > benchmarks/${{column}}.tsv 2>> {log}
        done
        """
