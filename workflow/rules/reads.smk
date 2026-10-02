rule sra_download:
    output:
        r1="resources/reads/{sample}_1.fastq.gz",
        r2="resources/reads/{sample}_2.fastq.gz",
    log:
        "logs/sra/{sample}.log",
    conda:
        "../envs/sra.yaml"
    threads: 4
    resources:
        runtime=120,
    params:
        accession=lambda wc: SAMPLES.loc[wc.sample, "sra"],
        outdir=lambda wc: f"resources/reads/{wc.sample}_tmp",
    shell:
        """
        (
            fasterq-dump {params.accession} --split-files --threads {threads} -O {params.outdir}
            gzip -c {params.outdir}/{params.accession}_1.fastq >{output.r1}
            gzip -c {params.outdir}/{params.accession}_2.fastq >{output.r2}
            rm -rf {params.outdir}
        ) >{log} 2>&1
        """


rule fastp:
    input:
        unpack(raw_reads),
    output:
        r1=temp("results/trimmed/{sample}_1.fastq.gz"),
        r2=temp("results/trimmed/{sample}_2.fastq.gz"),
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/fastp.yaml"
    threads: 4
    params:
        extra=config["fastp"]["extra"],
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--thread {threads} {params.extra} "
        "--html {output.html} --json {output.json} "
        "--report_title '{wildcards.sample}' > {log} 2>&1"
