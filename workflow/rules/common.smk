import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

SAMPLES = (
    pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    .fillna("")
    .set_index("sample", drop=False)
)
validate(SAMPLES, schema="../schemas/samples.schema.yaml")

CONTRASTS = list(config["diffexp"]["contrasts"])
REF = config["reference"]["name"]


wildcard_constraints:
    sample="|".join(SAMPLES.index),
    ref=REF.replace(".", r"\."),


def is_local(sample):
    return SAMPLES.loc[sample, "fq1"] != ""


def raw_reads(wildcards):
    """Local FASTQ files if given in the sample sheet, otherwise SRA downloads."""
    if is_local(wildcards.sample):
        return {
            "r1": SAMPLES.loc[wildcards.sample, "fq1"],
            "r2": SAMPLES.loc[wildcards.sample, "fq2"],
        }
    return {
        "r1": f"resources/reads/{wildcards.sample}_1.fastq.gz",
        "r2": f"resources/reads/{wildcards.sample}_2.fastq.gz",
    }


def ref_source(key):
    """URL references are downloaded into resources/, local paths are used as-is."""
    path = config["reference"][key]
    return path if not path.startswith(("http://", "https://", "ftp://")) else None


def ref_input(key, ext):
    local = ref_source(key)
    return local if local else f"resources/reference/{REF}.{ext}.download.gz"
