rule star_align:
    input:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz",
        index=f"resources/star_index/{REF}",
    output:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log_final="results/star/{sample}/Log.final.out",
    log:
        "logs/star/{sample}.log",
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        prefix="results/star/{sample}/",
        extra=config["star"]["extra"],
    shell:
        "STAR --runMode alignReads --runThreadN {threads} "
        "--genomeDir {input.index} --readFilesIn {input.r1} {input.r2} "
        "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate "
        "--outSAMattributes NH HI AS nM --outTmpDir {resources.tmpdir}/star_{wildcards.sample} "
        "--outFileNamePrefix {params.prefix} {params.extra} > {log} 2>&1"


rule samtools_index:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools/index_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} 2> {log}"


rule samtools_stats:
    input:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/qc/samtools/{sample}.stats.txt",
    log:
        "logs/samtools/stats_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input.bam} > {output} 2> {log}"
