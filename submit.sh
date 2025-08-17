#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2
#$ -j y
#$ -R y

WORKFLOW_NAME=$1
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
snakemake -s "$SCRIPT_DIR/workflows/$WORKFLOW_NAME" \
    --conda-prefix "/SAN/vyplab/alb_projects/" \
    --use-conda \
    --use-singularity \
    --singularity-args "--bind /SAN/vyplab/:/SAN/vyplab/,/scratch0/:/scratch0/" \
    --jobscript scripts/cluster_qsub.sh \
    --cluster-config scripts/cluster.yaml \
    --cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o submissions/{rule}.{wildcards}.\$JOB_ID.sh.o {cluster.submission_string}" \
    -j 8 \
    --rerun-triggers mtime \
    --nolock \
    --rerun-incomplete \
    --conda-frontend conda \
    --latency-wait 200
