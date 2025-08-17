#star - align to human genome

GENOME_DIR_HUMAN = config.get("genomeDir_human", config.get("dir_name_star_human", "STAR_index_human/"))
OUTFILEPREFIX_HUMAN = config.get("outFileNamePrefix_human", "human")

if IS_PAIRED:
    rule star_human:
        input:
            fq1 = "results/qc/{sample}_1.fq.gz",
            fq2 = "results/qc/{sample}_2.fq.gz"
        output:
            bam = "results/human/{sample}.Aligned.sortedByCoord.out.bam",
            log = "results/human/{sample}.Log.final.out"
        log:
            "logs/{sample}_star_human.log"
        params:
            prefix = lambda wildcards: f"results/human/{wildcards.sample}.",
            genome = GENOME_DIR_HUMAN
        threads: 4
        singularity:
            "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
        shell:
            """
            mkdir -p results/human
            STAR --runMode alignReads --genomeDir {params.genome} \
                 --readFilesIn {input.fq1} {input.fq2} \
                 --readFilesCommand zcat --runThreadN {threads} \
                 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
                 --outFileNamePrefix {params.prefix} 2> {log}
            """
else:
    rule star_human:
        input:
            fq = "results/qc/{sample}.fq.gz"
        output:
            bam = "results/human/{sample}.Aligned.sortedByCoord.out.bam",
            log = "results/human/{sample}.Log.final.out"
        log:
            "logs/{sample}_star_human.log"
        params:
            prefix = lambda wildcards: f"results/human/{wildcards.sample}.",
            genome = GENOME_DIR_HUMAN
        threads: 4
        singularity:
            "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
        shell:
            """
            mkdir -p results/human
            STAR --runMode alignReads --genomeDir {params.genome} \
                 --readFilesIn {input.fq} \
                 --readFilesCommand zcat --runThreadN {threads} \
                 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
                 --outFileNamePrefix {params.prefix} 2> {log}
            """