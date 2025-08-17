import re
import pandas as pd

configfile: "scripts/config.yaml"

df = pd.read_csv(config["samples_csv"], dtype=str)
sample_names = df["sample"].dropna().astype(str).tolist()
SAMPLE_RE = "|".join(map(re.escape, sample_names))

IS_PAIRED = config.get("end_type", "se") == "pe"
FILENAME_OUTPUT = config.get("filename_output", "VIRTUS.output.tsv")

def get_run_id(sample):
    return df.loc[df["sample"] == sample, "run"].values[0]

include: "../rules/download_sra.smk"
include: "../rules/fastp.smk"
include: "../rules/star_human.smk"
include: "../rules/extract_unmapped.smk"
include: "../rules/bam_to_fastq.smk"
include: "../rules/kz_filter.smk"
include: "../rules/star_virus.smk"
include: "../rules/filter_polyx.smk"
include: "../rules/samtools_coverage.smk"
include: "../rules/final_report.smk"

rule all:
    input:
        expand(f"results/final/{{sample}}_{FILENAME_OUTPUT}", sample=sample_names)
