#extract unmapped reads from star - human

wildcard_constraints:
    sample="|".join(map(re.escape, sample_names))

rule extract_unmapped_bam:
    input:
        bam = bam_for_sample
    output:
        bam = "results/unmapped/raw/{sample}_unmapped.bam"
    log:
        "logs/{sample}_extract_unmapped.log"
    threads: 4
    singularity:
        "docker://yyasumizu/bam_filter_polyx:1.3"
    shell:
        """
        mkdir -p results/unmapped/raw

        samtools view -h -@ {threads} -f 4 {input} | grep -v "uT:A:3" | \
        samtools view -@ {threads} -bS - > {output.bam} || touch {output.bam}

        samtools index {output.bam} || true
        """