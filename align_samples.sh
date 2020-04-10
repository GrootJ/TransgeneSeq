#!/bin/bash
#$ -cwd
#$ -l h_rt=1:00:00
#$ -pe node 8

#fed by align_samples_qsub.py
#if cluster is busy, time requirement can be lowered. time per job is usually ~5 minutes
#name of sample/job and qsub stdour/err are provided by align_samples_qsub.py
#sourcing bashrc seems necessary to load modules
source /etc/bashrc

module load java/jre/1.8.0_66
module load python/2.7.11
module load picard/2.6.0
module load bwa/0.7.15
module load samtools/1.3.1

python align_sample.py --fastq1 $fileR1 --fastq2 $fileR2 --ref_fasta $ref_fasta --out_dir $out_dir
