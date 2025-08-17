#calculate coverage using samtools

rule samtools_coverage:
    input:
        input_bam = "results/virus/{sample}.filtered.out.bam"
    output:
        output = "results/coverage/{sample}_virus.coverage.txt"
    log:
        "logs/{sample}.coverage.log"
    threads: 4
    singularity:
        "docker://quay.io/biocontainers/samtools:1.15--h1170115_1"
    shell:
        """
        mkdir results/coverage
        samtools coverage -@ {threads} {input.input_bam} > {output.output} 2> {log}
        """
