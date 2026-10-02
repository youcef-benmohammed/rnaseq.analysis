rule multiqc:
    input:
        expand("results/qc/fastp/{sample}.json", sample=SAMPLES.index),
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES.index),
        expand("results/qc/samtools/{sample}.stats.txt", sample=SAMPLES.index),
        "results/counts/featureCounts.txt.summary",
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --outdir results/multiqc --filename multiqc_report.html "
        "{input} > {log} 2>&1"
