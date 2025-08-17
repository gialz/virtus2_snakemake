# KZ filter - filter low complexity

wildcard_constraints:
    sample="|".join(map(re.escape, sample_names))

KZ_THRESHOLD = config.get("kz_threshold", 0.1)

if IS_PAIRED:
    rule kz_filter_pe:
        wildcard_constraints:
            sample="(?!.*_kz$)(?!.*_unmapped$)[^/]+"
        input:
            fq1 = "results/unmapped/{sample}_1.fq",
            fq2 = "results/unmapped/{sample}_2.fq"
        output:
            fq1_out = "results/unmapped/{sample}_kz_1.fq",
            fq2_out = "results/unmapped/{sample}_kz_2.fq"
        log:
            "logs/{sample}_kz.log"
        params:
            threshold = KZ_THRESHOLD
        singularity:
            "docker://yyasumizu/ko:0.1"
        shell:
            """
            mkdir -p results/unmapped logs
            kz --filter --threshold {params.threshold} < {input.fq1} > {output.fq1_out} 2> {log}
            kz --filter --threshold {params.threshold} < {input.fq2} > {output.fq2_out} 2>> {log}
            """

    rule fastq_pair:
        wildcard_constraints:
            sample="(?!.*_kz$)(?!.*_unmapped$)[^/]+"
        input:
            fq1 = "results/unmapped/{sample}_kz_1.fq",
            fq2 = "results/unmapped/{sample}_kz_2.fq"
        output:
            fq1_paired = "results/unmapped/{sample}_kz_1.fq.paired.fq",
            fq2_paired = "results/unmapped/{sample}_kz_2.fq.paired.fq"
        log:
            "logs/{sample}_fq.paired.log"
        singularity:
            "docker://quay.io/biocontainers/fastq-pair:1.0--he1b5a44_1"
        shell:
            """
            mkdir -p results/unmapped logs
            fastq_pair {input.fq1} {input.fq2} 2> {log}
            """

else:
    rule kz_filter_se:
        wildcard_constraints:
            sample="(?!.*_kz$)(?!.*_unmapped$)[^/]+"
        input:
            fq = "results/unmapped/{sample}.fq"
        output:
            fq_out = "results/unmapped/{sample}_kz.fq"
        log:
            "logs/{sample}_kz.log"
        params:
            threshold = KZ_THRESHOLD
        singularity:
            "docker://yyasumizu/ko:0.1"
        shell:
            """
            kz --filter --threshold {params.threshold} < {input.fq} > {output.fq_out} 2> {log}
            """