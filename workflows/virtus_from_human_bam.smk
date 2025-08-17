import pandas as pd

configfile: "scripts/config.yaml"

bam_df = pd.read_csv(config["samples_csv"]).dropna(subset=["sample", "bam"])
bam_df["sample"] = bam_df["sample"].astype(str)
bam_df["bam"] = bam_df["bam"].astype(str)
bam_map = dict(zip(bam_df["sample"], bam_df["bam"]))
sample_names = list(bam_map.keys()) 

IS_PAIRED = config.get("end_type", "se") == "pe"
FILENAME_OUTPUT = config.get("filename_output", "VIRTUS.output.tsv")

def bam_for_sample(wc):
    bam_path = bam_map.get(wc.sample, "").strip()
    if bam_path and os.path.exists(bam_path):
        return bam_path
    # fallback to default
    return f"results/unmapped/raw/{wc.sample}_unmapped.bam"

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
