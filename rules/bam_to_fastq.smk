#bam to fastq

wildcard_constraints:
    sample="|".join(map(re.escape, sample_names))

if IS_PAIRED:
    rule bam_to_fastq_pe:
        wildcard_constraints:
            sample="(?!.*_kz$)(?!.*_unmapped$)[^/]+"
        input:
            bam = "results/unmapped/raw/{sample}_unmapped.bam"
        output:
            fq1 = "results/unmapped/{sample}_1.fq",
            fq2 = "results/unmapped/{sample}_2.fq"
        log:
            "logs/{sample}_bamtofastq_pe.log"
        threads: 2
        singularity:
            "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
        shell:
            """
            bedtools bamtofastq -i {input.bam} \
                -fq {output.fq1} -fq2 {output.fq2} || touch {output.fq1} {output.fq2}
            """

else:
    rule bam_to_fastq_se:
        wildcard_constraints:
            sample="(?!.*_kz$)(?!.*_unmapped$)[^/]+"
        input:
            bam = "results/unmapped/raw/{sample}_unmapped.bam"
        output:
            fq_single = "results/unmapped/{sample}.fq"
        log:
            "logs/{sample}_bamtofastq_se.log"
        threads: 2
        singularity:
            "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
        shell:
            """
            bedtools bamtofastq -i {input.bam} \
                -fq {output.fq_single} || touch {output.fq_single}
           """