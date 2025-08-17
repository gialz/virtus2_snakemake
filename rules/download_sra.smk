# sra download

if IS_PAIRED:
    rule download_sra_pe:
        output:
            r1 = "data/{sample}_1.fastq.gz",
            r2 = "data/{sample}_2.fastq.gz"
        log:
            "logs/{sample}_download.log"
        params:
            sra_toolkit_path = config.get("sra_toolkit", ""),
            run_id = lambda wildcards: get_run_id(wildcards.sample)
        threads: 4
        shell:
            """
            mkdir -p data logs
            if [ -f {output.r1} ] && [ -f {output.r2} ]; then
                echo "FASTQs exist for {params.run_id}, skipping." >> {log}
            else
                {params.sra_toolkit_path}prefetch --max-size 70000000000 {params.run_id} 2>> {log}
                {params.sra_toolkit_path}fasterq-dump -e {threads} -x {params.run_id} --split-files -O data 2>> {log}
                gzip -f data/{params.run_id}_1.fastq data/{params.run_id}_2.fastq 2>> {log}
                rm -rf {params.run_id}
            fi
            """

else:
    rule download_sra_se:
        output:
            fq = "data/{sample}.fastq.gz"
        log:
            "logs/{sample}_download.log"
        params:
            sra_toolkit_path = config.get("sra_toolkit", ""),
            run_id = lambda wildcards: get_run_id(wildcards.sample)
        threads: 4
        shell:
            """
            mkdir -p data logs
            if [ -f {output.fq} ]; then
                echo "FASTQ exists for {params.run_id}, skipping." >> {log}
            else
                {params.sra_toolkit_path}prefetch --max-size 70000000000 {params.run_id} 2>> {log}
                {params.sra_toolkit_path}fasterq-dump -e {threads} -x {params.run_id} -O data 2>> {log}
                gzip -f data/{params.run_id}.fastq 2>> {log}
                rm -rf {params.run_id}
            fi
            """

