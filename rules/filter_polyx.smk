# filter poly-x tails

rule bam_filter_polyx:
    input:
        input_bam = "results/virus/{sample}_virus.Aligned.sortedByCoord.out.bam"
    output:
        output = "results/virus/{sample}.filtered.out.bam"
    threads: 4
    singularity:
        "docker://yyasumizu/bam_filter_polyx:1.3"
    shell:
        """
        mkdir -p results/virus
        samtools view -@ {threads} -h {input.input_bam} \
        | grep -Fv "AAAAAAAAAAAAAAAAAAAA" \
        | grep -Fv "TTTTTTTTTTTTTTTTTTTT" \
        | grep -Fv "TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG" \
        | samtools view -@ {threads} -b > {output.output}
        samtools index -@ {threads} {output.output}
        """
