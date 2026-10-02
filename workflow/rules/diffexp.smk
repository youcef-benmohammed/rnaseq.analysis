rule deseq2:
    input:
        counts="results/counts/counts.tsv",
        samples=config["samples"],
    output:
        results=expand("results/diffexp/{contrast}.results.tsv", contrast=CONTRASTS),
        volcano=expand("results/diffexp/{contrast}.volcano.png", contrast=CONTRASTS),
        ma=expand("results/diffexp/{contrast}.ma.png", contrast=CONTRASTS),
        normalized="results/diffexp/normalized_counts.tsv",
        pca="results/diffexp/pca.png",
        heatmap="results/diffexp/sample_distances.png",
        rds="results/diffexp/dds.rds",
    log:
        "logs/deseq2.log",
    conda:
        "../envs/deseq2.yaml"
    threads: 2
    params:
        diffexp=config["diffexp"],
        contrasts=CONTRASTS,
    script:
        "../scripts/deseq2.R"
