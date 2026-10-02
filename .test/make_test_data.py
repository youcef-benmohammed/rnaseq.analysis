"""Generate a tiny synthetic RNA-seq dataset for CI and quick tests.

Genome: 2 chromosomes, 40 multi-exon genes (GTF). Paired-end 2x100 bp reads
are sampled from transcripts for 2 conditions x 3 replicates; genes 1-8 are
4x up and genes 9-16 are 4x down in condition TREAT, so the full workflow
(STAR -> featureCounts -> DESeq2) has a known answer.

    python .test/make_test_data.py   # writes .test/data/ and .test/config/samples.tsv
    snakemake --directory .test --cores 4 --sdm conda   # run the workflow on it
"""
import gzip
import random
from pathlib import Path

random.seed(42)
OUT = Path(__file__).parent / "data"
OUT.mkdir(parents=True, exist_ok=True)

N_PAIRS = 3000
READ_LEN = 100
FRAG_MEAN, FRAG_SD = 250, 25
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


# --- genome and annotation ----------------------------------------------------
chroms = {f"chr{i}": "".join(random.choice("ACGT") for _ in range(60_000)) for i in (1, 2)}
genes = []  # (gene_id, chrom, strand, [(start, end), ...]) 1-based inclusive
gid = 0
for chrom, seq in chroms.items():
    pos = 1_000
    while pos < len(seq) - 3_000 and gid < 40:
        gid += 1
        exons, cur = [], pos
        for _ in range(random.randint(1, 3)):
            length = random.randint(300, 700)
            exons.append((cur, cur + length - 1))
            cur += length + random.randint(150, 400)  # intron
        genes.append((f"GENE{gid:03d}", chrom, random.choice("+-"), exons))
        pos = cur + random.randint(300, 600)

with open(OUT / "genome.fa", "w") as fa:
    for chrom, seq in chroms.items():
        fa.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            fa.write(seq[i:i + 60] + "\n")

with open(OUT / "annotation.gtf", "w") as gtf:
    for g, chrom, strand, exons in genes:
        attrs = f'gene_id "{g}"; transcript_id "{g}.1"; gene_biotype "protein_coding";'
        gtf.write(f"{chrom}\tsynthetic\tgene\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n")
        gtf.write(f"{chrom}\tsynthetic\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n")
        for n, (s, e) in enumerate(exons, 1):
            gtf.write(f"{chrom}\tsynthetic\texon\t{s}\t{e}\t.\t{strand}\t.\t{attrs} exon_number \"{n}\";\n")

transcripts = {g: "".join(chroms[c][s - 1:e] for s, e in ex) for g, c, _, ex in genes}

# --- expression and reads -------------------------------------------------------
base = {g: random.lognormvariate(0, 0.8) for g in transcripts}
fold = {g: 4.0 if i < 8 else 0.25 if i < 16 else 1.0 for i, g in enumerate(transcripts)}


def simulate(sample, condition):
    weights = []
    for g, t in transcripts.items():
        w = base[g] * len(t) * (fold[g] if condition == "TREAT" else 1.0)
        weights.append(w * random.lognormvariate(0, 0.15))  # biological noise
    names = list(transcripts)
    with gzip.open(OUT / f"{sample}_R1.fastq.gz", "wt") as r1, \
         gzip.open(OUT / f"{sample}_R2.fastq.gz", "wt") as r2:
        for i in range(N_PAIRS):
            t = transcripts[random.choices(names, weights)[0]]
            frag = max(READ_LEN, min(len(t), int(random.gauss(FRAG_MEAN, FRAG_SD))))
            start = random.randint(0, len(t) - frag)
            f = t[start:start + frag]
            if random.random() < 0.5:  # unstranded library
                f = revcomp(f)
            qual = "I" * READ_LEN
            r1.write(f"@{sample}.{i}/1\n{f[:READ_LEN]}\n+\n{qual}\n")
            r2.write(f"@{sample}.{i}/2\n{revcomp(f)[:READ_LEN]}\n+\n{qual}\n")


rows = ["sample\tcondition\treplicate\tsra\tfq1\tfq2"]
for condition in ("CTRL", "TREAT"):
    for rep in (1, 2, 3):
        s = f"{condition}_rep{rep}"
        simulate(s, condition)
        rows.append(f"{s}\t{condition}\t{rep}\t\tdata/{s}_R1.fastq.gz\tdata/{s}_R2.fastq.gz")
(Path(__file__).parent / "config" / "samples.tsv").write_text("\n".join(rows) + "\n")
print(f"{len(genes)} genes, {N_PAIRS} read pairs x 6 samples written to {OUT}")
