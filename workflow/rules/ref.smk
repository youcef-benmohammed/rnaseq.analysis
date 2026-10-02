rule download_reference:
    output:
        temp("resources/reference/{ref}.{ext}.download.gz"),
    log:
        "logs/reference/download_{ref}.{ext}.log",
    wildcard_constraints:
        ext="fa|gtf",
    conda:
        "../envs/utils.yaml"
    params:
        url=lambda wc: config["reference"]["fasta" if wc.ext == "fa" else "gtf"],
    shell:
        "curl -fsSL --retry 3 -o {output} {params.url} 2> {log}"


rule prepare_reference:
    """Decompress (if needed) so that STAR and featureCounts get plain text."""
    input:
        lambda wc: ref_input("fasta" if wc.ext == "fa" else "gtf", wc.ext),
    output:
        "resources/reference/{ref}.{ext}",
    log:
        "logs/reference/prepare_{ref}.{ext}.log",
    wildcard_constraints:
        ext="fa|gtf",
    conda:
        "../envs/utils.yaml"
    shell:
        "(if gzip -t {input} 2>/dev/null; then gzip -dc {input}; else cat {input}; fi) > {output} 2> {log}"


rule star_index:
    input:
        fasta=f"resources/reference/{REF}.fa",
        gtf=f"resources/reference/{REF}.gtf",
    output:
        directory(f"resources/star_index/{REF}"),
    log:
        f"logs/star/index_{REF}.log",
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        sa=config["star"]["sa_index_nbases"],
        tmp=lambda wc, resources: f"{resources.tmpdir}/star_index_{REF}",
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate --runThreadN {threads} "
        "--genomeDir {output} --genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} --genomeSAindexNbases {params.sa} "
        "--outTmpDir {params.tmp} "
        "--outFileNamePrefix {output}/ > {log} 2>&1"
