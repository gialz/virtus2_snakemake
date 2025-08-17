# fastp - quality control

FASTP_LENGTH = config.get("fastp_length", 40)

if IS_PAIRED:
    rule fastp_pe:
        input:
            fq1 = lambda wildcards: get_fastqs(wildcards.sample)[0],
            fq2 = lambda wildcards: get_fastqs(wildcards.sample)[1]
        output:
            fq1_out = "results/qc/{sample}_1.fq.gz",
            fq2_out = "results/qc/{sample}_2.fq.gz",
            html    = "results/qc/{sample}_fastp_report.html",
            json    = "results/qc/{sample}_fastp_report.json"
        log:
            "logs/{sample}_fastp.log"
        params:
            length = FASTP_LENGTH
        threads: 8
        singularity:
            "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
        shell:
            """
            mkdir -p results/qc logs
            fastp --in1 {input.fq1} --in2 {input.fq2} \
                  --out1 {output.fq1_out} --out2 {output.fq2_out} \
                  --thread {threads} --trim_poly_x \
                  --html {output.html} --json {output.json} \
                  --length_required {params.length} 2> {log}
            """

else:
    rule fastp_se:
        input:
            fq = lambda wildcards: get_fastqs(wildcards.sample)[0]
        output:
            fq_out = "results/qc/{sample}.fq.gz",
            html   = "results/qc/{sample}_fastp_report.html",
            json   = "results/qc/{sample}_fastp_report.json"
        log:
            "logs/{sample}_fastp.log"
        params:
            length = FASTP_LENGTH
        threads: 8
        singularity:
            "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
        shell:
            """
            mkdir -p results/qc logs
            fastp --in1 {input.fq} --out1 {output.fq_out} \
                  --thread {threads} --trim_poly_x \
                  --html {output.html} --json {output.json} \
                  --length_required {params.length} 2> {log}
            """

