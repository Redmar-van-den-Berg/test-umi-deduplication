include: "common.smk"


rule all:
    input:
        concat=expand(
            "humid/{sample}/forward_dedup.fastq.gz",
            sample=samples,
        ),
        bam=expand("{sample}/umi-tools/{sample}.bam", sample=samples),
        stats="umi-stats.tsv",
        fastq=expand("{sample}/umi-tools/forward.fastq.gz", sample=samples),
        humid_after_umi_tools=expand(
            "{sample}/umi-tools/humid/forward_dedup.fastq.gz", sample=samples
        ),
        umi_tools_after_humid=expand(
            "{sample}/humid/umi-tools/{sample}.bam", sample=samples
        ),
        multiqc="multiqc_report.html",
        multiqc_humid_after_umi_tools="multiqc_report_humid_after_umi_tools.html",
        mutliqc_umi_tools_after_humid="multiqc_report_umi_tools_after_humid.html",
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


rule fastqc:
    """Runs FastQC"""
    input:
        fq1=rules.concat.output.forw,
        fq2=rules.concat.output.rev,
    output:
        folder=directory("{sample}/concat/fastqc_{sample}"),
        forw="{sample}/concat/fastqc_{sample}/forward_fastqc.zip",
        rev="{sample}/concat/fastqc_{sample}/reverse_fastqc.zip",
    params:
        xms="4096M",
        xmx="4096M",
        fastqc_dir="usr/local/opt/fastqc-0.11.9",
    log:
        "log/fastqc_raw.{sample}.txt",
    threads: 4
    container:
        containers["fastqc"]
    shell:
        """
        mkdir -p {output.folder} tmp

        FASTQC_DIR=/{params.fastqc_dir}
        export CLASSPATH="$FASTQC_DIR:$FASTQC_DIR/sam-1.103.jar:$FASTQC_DIR/jbzip2-0.9.jar:$FASTQC_DIR/cisd-jhdf5.jar"

        java -Djava.awt.headless=true -Xms{params.xms} -Xmx{params.xmx} \
            -Dfastqc.output_dir={output.folder} \
            -Dfastqc.io.tmpdir=tmp \
            -Dfastqc.unzip=true \
            -Dfastqc.nogroup=true \
            -Dfastqc.threads={threads} \
            uk.ac.babraham.FastQC.FastQCApplication \
            {input.fq1:q} {input.fq2:q} 2> {log}
        """


rule humid:
    """Run HUMID on the fastq files"""
    input:
        forw=rules.concat.output.forw,
        rev=rules.concat.output.rev,
        umi=rules.concat.output.umi,
    params:
        cluster_method="-e" if config["cluster_method"] == "maximum" else "",
        stack_size_kb=102400,
    output:
        forw="humid/{sample}/forward_dedup.fastq.gz",
        rev="humid/{sample}/reverse_dedup.fastq.gz",
        umi="humid/{sample}/umi_dedup.fastq.gz",
        stats="humid/{sample}/stats.dat",
        clusters="humid/{sample}/clusters.dat",
        counts="humid/{sample}/counts.dat",
        neigh="humid/{sample}/neigh.dat",
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

        ulimit -Ss {params.stack_size_kb}

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
        main="{sample}/align/Log.out",
        progress="{sample}/align/Log.progress.out",
        final="{sample}/align/Log.final.out",
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
        bam="{sample}/umi-tools/{sample}.bam",
    log:
        "log/{sample}_umi_tools.log",
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
        umi_tools=expand("log/{sample}_umi_tools.log", sample=samples),
        humid=expand("humid/{sample}/stats.dat", sample=samples),
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
                --humid {input.humid} \
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


# Run HUMID on the output of umi-tools
use rule humid as humid_after_umi_tools with:
    input:
        forw=rules.bam_to_fastq.output.forw,
        rev=rules.bam_to_fastq.output.rev,
        umi=rules.bam_to_fastq.output.umi,
    output:
        forw="{sample}/umi-tools/humid/forward_dedup.fastq.gz",
        rev="{sample}/umi-tools/humid/reverse_dedup.fastq.gz",
        umi="{sample}/umi-tools/humid/umi_dedup.fastq.gz",
        stats="{sample}/umi-tools/humid/stats.dat",
    log:
        "log/{sample}.humid_after_umi_tools.txt",
    benchmark:
        repeat("benchmarks/humid_after_umi_tools_{sample}.tsv", config["repeat"])


# Add the UMI's to the read name after running HUMID
use rule add_umi as add_umi_after_humid with:
    input:
        forw=rules.humid.output.forw,
        rev=rules.humid.output.rev,
        umi=rules.humid.output.umi,
        add_umi=srcdir("scripts/add-umi.py"),
    output:
        forw="{sample}/humid/umi-tools/{sample}.umi.forward.fastq.gz",
        rev="{sample}/humid/umi-tools/{sample}.umi.reverse.fastq.gz",
    log:
        "log/{sample}.add_umi_after_humid.txt",


# Align the reads to the reference
use rule align_vars as align_after_humid with:
    input:
        fq1=rules.add_umi_after_humid.output.forw,
        fq2=rules.add_umi_after_humid.output.rev,
        index=config["star_index"],
        gtf=config["gtf"],
    output:
        bam="{sample}/humid/umi-tools/align/Aligned.sortedByCoord.out.bam",
    log:
        "log/{sample}.align_after_humid.txt",
    benchmark:
        repeat("benchmarks/STAR_after_humid_{sample}.tsv", config["repeat"])


use rule index_bamfile as index_after_humid with:
    input:
        bam=rules.align_after_humid.output.bam,
    output:
        bai="{sample}/humid/umi-tools/align/Aligned.sortedByCoord.out.bam.bai",
    log:
        "log/{sample}.index_after_humid.txt",


# Run umi-tools
use rule umi_tools as umi_tools_after_humid with:
    input:
        bam=rules.align_after_humid.output.bam,
        bai=rules.index_after_humid.output.bai,
    output:
        bam="{sample}/humid/umi-tools/{sample}.bam",
    log:
        "log/{sample}_umi_tools_after_humid.log",
    benchmark:
        repeat("benchmarks/umi_tools_after_humid_{sample}.tsv", config["repeat"])


# Run MultiQC on HUMID output
rule multiqc:
    input:
        stats=get_log_files(),
        config=srcdir("cfg/multiqc.yml"),
    params:
        filelist="multiqc_filelist.txt",
        depth=2,
    output:
        html="multiqc_report.html",
    log:
        "log/multiqc.txt",
    container:
        containers["multiqc"]
    shell:
        """
        rm -f {params.filelist}

        for fname in {input.stats}; do
            echo $fname >> {params.filelist}
        done

        multiqc \
        --force \
        --dirs \
        --dirs-depth {params.depth} \
        --fullnames \
        --fn_as_s_name \
        --file-list {params.filelist} \
        --config {input.config} \
        --filename {output.html} 2> {log}
        """


# Run MultiQC on HUMID output after deduplication with UMI-Tools
use rule multiqc as multiqc_humid_after_umi_tools with:
    input:
        stats=get_humid_after_umi_tools(),
        config=srcdir("cfg/multiqc_humid_after_umi_tools.yml"),
    params:
        filelist="multiqc_humid_after_umi_tools_filelist.txt",
        depth=3,
    output:
        html="multiqc_report_humid_after_umi_tools.html",
    log:
        "log/multiqc_humid_after_umi_tools.txt",


# Run MultiQC on UMI-Tools output after deduplication with HUMID
use rule multiqc as multiqc_umi_tools_after_humid with:
    input:
        stats=get_umi_tools_after_humid(),
        config=srcdir("cfg/multiqc_umi_tools_after_humid.yml"),
    params:
        filelist="multiqc_umi_tools_after_humid_filelist.txt",
        depth=1,
    output:
        html="multiqc_report_umi_tools_after_humid.html",
    log:
        "log/multiqc_umi_tools_after_humid.txt",


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


rule picard_multiple_stats:
    input:
        bam=rules.align_vars.output.bam,
        bai=rules.index_bamfile.output.bai,
        ref=config["reference"],
        ref_dict=config["reference_dict"],
    params:
        prefix=lambda wildcards, output: output.alignment_summary.split(".")[0],
    output:
        alignment_summary="{sample}/align/multiple_metrics.alignment_summary_metrics",
        inserts="{sample}/align/multiple_metrics.insert_size_metrics",
        base_dist="{sample}/align/multiple_metrics.base_distribution_by_cycle_metrics",
        qual="{sample}/align/multiple_metrics.quality_by_cycle_metrics",
        qual_dist="{sample}/align/multiple_metrics.quality_distribution_metrics",
    log:
        "log/picard_multiple_stats.{sample}.txt",
    threads: 1
    container:
        containers["picard"]
    shell:
        """
        picard -Xmx4G  CollectMultipleMetrics \
            VALIDATION_STRINGENCY=LENIENT \
            R={input.ref} \
            I={input.bam} \
            O={params.prefix} 2> {log}
        """
