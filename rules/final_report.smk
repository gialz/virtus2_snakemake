#final report
import os

FILENAME_OUTPUT = config.get("filename_output", "VIRTUS.output.tsv")

rule final_report:
    input:
        star_log = bam_for_sample,
        coverage = "results/coverage/{sample}_virus.coverage.txt"
    output:
        "results/final/{sample}_" + FILENAME_OUTPUT
    run:
        import pandas as pd
        import os

        # Read STAR log safely
        df_log = pd.read_csv(input.star_log, sep='\t', header=None, index_col=0)
        df_log.index = df_log.index.str.strip()  # normalize keys
        num_reads = int(df_log.loc['Uniquely mapped reads number', 1]) + \
                    int(df_log.loc['Number of reads mapped to multiple loci', 1])

        # Read coverage
        df_cov = pd.read_csv(input.coverage, sep='\t')
        df_cov.columns = [
            'virus', 'startpos', 'endpos', 'numreads', 'covbases',
            'coverage', 'meandepth', 'meanbaseq', 'meanmapq'
        ]

        # Adjust for paired-end
        if IS_PAIRED:
            df_cov['numreads'] /= 2

        df_cov['rate_hit'] = df_cov['numreads'] / num_reads
        df_cov = df_cov[df_cov['numreads'] > 0]
        df_cov = df_cov.sort_values(by='rate_hit', ascending=False)

        os.makedirs("results/final", exist_ok=True)
        df_cov.to_csv(output[0], index=False, sep='\t')
