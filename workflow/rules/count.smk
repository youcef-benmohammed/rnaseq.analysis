rule featurecounts:
    """Gene-level counts. Duplicates are NOT removed: in RNA-seq most duplicates
    come from highly expressed genes, and removing them biases quantification."""
    input:
        bams=expand(
            "results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES.index
        ),
        gtf=f"resources/reference/{REF}.gtf",
    output:
        counts="results/counts/featureCounts.txt",
        summary="results/counts/featureCounts.txt.summary",
    log:
        "logs/featurecounts.log",
    conda:
        "../envs/subread.yaml"
    threads: 4
    params:
        strand=config["featurecounts"]["strandedness"],
        extra=config["featurecounts"]["extra"],
    shell:
        "featureCounts -p --countReadPairs -B -C -s {params.strand} -T {threads} "
        "-a {input.gtf} -o {output.counts} {params.extra} {input.bams} > {log} 2>&1"


rule count_matrix:
    """Clean gene x sample matrix (sample names instead of BAM paths)."""
    input:
        "results/counts/featureCounts.txt",
    output:
        "results/counts/counts.tsv",
    log:
        "logs/count_matrix.log",
    conda:
        "../envs/python.yaml"
    params:
        samples=list(SAMPLES.index),
    script:
        "../scripts/count_matrix.py"
