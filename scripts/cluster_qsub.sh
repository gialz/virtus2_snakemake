#$ -S /bin/bash

#$ -cwd
#$ -V
# join stdout and stderr output
#$ -j y
#$ -sync y


echo $JOB_ID

{exec_job}