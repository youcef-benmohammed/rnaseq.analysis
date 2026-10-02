"""Turn featureCounts output into a clean gene x sample count matrix."""
import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

fc = pd.read_csv(snakemake.input[0], sep="\t", comment="#", index_col=0)
counts = fc.iloc[:, 5:]  # drop Chr, Start, End, Strand, Length
if counts.shape[1] != len(snakemake.params.samples):
    raise ValueError("Number of count columns does not match the number of samples")
counts.columns = snakemake.params.samples  # BAMs are passed in sample-sheet order
counts.index.name = "gene_id"
counts.to_csv(snakemake.output[0], sep="\t")
print(f"{counts.shape[0]} genes x {counts.shape[1]} samples", file=sys.stderr)
