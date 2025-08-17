import pandas as pd
import os

configfile: "scripts/config.yaml"

df = pd.read_csv(config["samples_csv"], dtype=str)

# Keep only rows where at least one FASTQ path is given
if config.get("end_type", "se") == "pe":
    df = df[df["fastq_1"].notna() & df["fastq_2"].notna()]
else:
    df = df[df["fastq_1"].notna()]

sample_names = df["sample"].astype(str).tolist()

IS_PAIRED = config.get("end_type", "se") == "pe"
FILENAME_OUTPUT = config.get("filename_output", "VIRTUS.output.tsv")

def get_fastqs(sample):
    row = df[df["sample"] == str(sample)].iloc[0]
    if IS_PAIRED:
        return [
            row["fastq_1"] if row["fastq_1"].strip() else f"data/{sample}_1.fastq.gz",
            row["fastq_2"] if row["fastq_2"].strip() else f"data/{sample}_2.fastq.gz"
        ]
    else:
        return [
            row["fastq_1"] if row["fastq_1"].strip() else f"data/{sample}.fastq.gz"
        ]

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
