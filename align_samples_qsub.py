#!/usr/bin/python
import sys
import re
import os
import argparse
import fnmatch
import subprocess

#loops through all fastqs in directory specified. file extension should be changed if need be, for gzipped files add .gz
#usage:
#module add python/2.7.11  #(module add python should work, as this is the default version)
#python align_samples_qsub.py --fastq_dir /../project/fastq --ref_fasta /../project/reference.fa --out_dir /../project --stdout_dir /../project/stdout 
# modify file extensions below to match fastq file names (e.g. fastq.gz instead of .fastq.gz)
#directories are made below if they do not exist already
def parseArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq_dir', required=True)
    parser.add_argument('--ref_fasta', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--stdout_dir', required=True)

    args = vars(parser.parse_args())

    return args


if __name__ == '__main__':
	
    args = parseArgs()
    fastq_dir = args['fastq_dir']
    ref_fasta = args['ref_fasta']
    out_dir = args['out_dir']
    qsub_stdout = args['stdout_dir']


    cmd = 'module add bwa/0.7.15; bwa index ' + ref_fasta
    subprocess.check_call(cmd, shell=True)

    #out_dir = os.path.join(out_dir, 'bam')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(qsub_stdout):
        os.makedirs(qsub_stdout)

    file_R1 = []
    file_R2 = []
    for file in os.listdir(fastq_dir):
        if fnmatch.fnmatch(file, '*R1_001.fastq.gz'):
            file_R1.append(os.path.join(fastq_dir, file))
        elif fnmatch.fnmatch(file, '*R2_001.fastq.gz'):
            file_R2.append(os.path.join(fastq_dir, file))
    file_R1.sort()
    file_R2.sort()

    if len(file_R1) != len(file_R2):
        print("Warnings: Not equal number of read1 and read2 fastq files!\n")
   
    for i in range(len(file_R1)):
        print(i)
        sample_name = re.sub('.R1_001.fastq.gz', '', os.path.basename(file_R1[i]))
        sample_name = 'Transgen' + sample_name
        cmd = 'qsub -v fileR1=' +  file_R1[i] + ' -v fileR2=' + file_R2[i] + ' -N ' + sample_name +  ' -v ref_fasta=' + ref_fasta + ' -v out_dir=' + out_dir + ' -o ' + qsub_stdout +' -e ' + qsub_stdout  + ' align_samples.sh'
        print(cmd)
        subprocess.check_call(cmd, shell=True)

    cmd = 'qsub -N wait_pro -hold_jid \"Transgen*\" -sync y -b y echo \"DONE!\"'
    print(cmd)
    subprocess.check_call(cmd, shell=True)


