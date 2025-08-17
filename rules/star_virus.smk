#star - align to virus genome

GENOME_DIR_VIRUS = config.get("genomeDir_virus", config.get("dir_name_star_virus", "STAR_index_virus"))

if IS_PAIRED:
    rule star_virus_pe:
        input:
            fq1 = "results/unmapped/{sample}_kz_1.fq.paired.fq",
            fq2 = "results/unmapped/{sample}_kz_2.fq.paired.fq"
        output:
            bam = "results/virus/{sample}_virus.Aligned.sortedByCoord.out.bam",
            log = "results/virus/{sample}_virus.Log.final.out"
        log:
            "logs/{sample}_star_virus.log"
        params:
            prefix = lambda wc: f"results/virus/{wc.sample}_virus.",
            genome = GENOME_DIR_VIRUS
        threads: 4
        singularity:
            "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
        shell:
            r"""
            mkdir -p results/virus logs
            STAR --runMode alignReads --genomeDir {params.genome} \
                 --readFilesIn {input.fq1} {input.fq2} \
                 --runThreadN {threads} --outSAMtype BAM SortedByCoordinate \
                 --outFileNamePrefix {params.prefix} 2> {log}
            """
else:
    rule star_virus_se:
        input:
            fq = "results/unmapped/{sample}_kz.fq"
        output:
            bam = "results/virus/{sample}_virus.Aligned.sortedByCoord.out.bam",
            log = "results/virus/{sample}_virus.Log.final.out"
        log:
            "logs/{sample}_star_virus.log"
        params:
            prefix = lambda wc: f"results/virus/{wc.sample}_virus.",
            genome = GENOME_DIR_VIRUS
        threads: 4
        singularity:
            "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
        shell:
            r"""
            mkdir -p results/virus logs
            STAR --runMode alignReads --genomeDir {params.genome} \
                 --readFilesIn {input.fq} \
                 --runThreadN {threads} --outSAMtype BAM SortedByCoordinate \
                 --outFileNamePrefix {params.prefix} 2> {log}
            """
