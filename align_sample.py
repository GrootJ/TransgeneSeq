#!/usr/bin/python
import sys
import re
import os
import subprocess
import argparse

#module load bwa/0.7.12

#performs alignment via bwa mem and soft clipping, sorting and marking of duplicates and read groups via picard
#usage: via align_samples_qsub.py
def parseArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq1', required=True)
    parser.add_argument('--fastq2', required=True)
    parser.add_argument('--ref_fasta', required=True)
    parser.add_argument('--out_dir', default='.')

    args = vars(parser.parse_args())

    return args
 

if __name__ == '__main__':

    args = parseArgs()

    sample_name = re.sub('.L001_R1_001.fastq.gz', '', os.path.basename(args['fastq1']))

    # Align
    print("Starting align...")
    sam_out = os.path.join(args['out_dir'], sample_name + '.sam')
    cmd = 'bwa mem ' + '-t 4 ' + '-O 6 ' + '-a ' + '-M ' + args['ref_fasta'] + ' ' + args['fastq1'] + ' ' + args['fastq2'] + ' > ' + sam_out
    subprocess.check_call(cmd, shell=True)

    # Clean, soft clip beyond of reference alignment and set unmapped read map to 0
    print("Starting clean...")
    clean_out = re.sub('.sam', '.clean.sam', sam_out)
    cmd = 'picard CleanSam INPUT=' + sam_out + ' OUTPUT=' + clean_out
    subprocess.check_call(cmd, shell=True)

    # Sort
    print("Starting sort...")
    sort_out = re.sub('.clean.sam', '.sorted.sam', clean_out)
    cmd = 'picard SortSam INPUT=' + clean_out + ' OUTPUT=' + sort_out + ' SORT_ORDER=coordinate CREATE_INDEX=false'    
    subprocess.check_call(cmd, shell=True)

    # Add read group ID
    print("Starting read group ID update...")
    rg_out = re.sub('.sorted.sam', '.sorted.rg.bam', sort_out)
    cmd = 'picard AddOrReplaceReadGroups INPUT=' + sort_out + ' OUTPUT=' + rg_out + ' SORT_ORDER=coordinate RGID=' + sample_name + ' RGLB=' + sample_name + ' RGPL=Illumina RGSM=' + sample_name + ' RGPU=' + sample_name + ' CREATE_INDEX=false' 
    subprocess.check_call(cmd, shell=True)

    # Mark Duplicates
    print("Marking Duplicates...")
    mark_dup_out = re.sub('.rg.bam', '.rg.dups_marked.bam', rg_out)
    metrics_out = re.sub('.rg.bam', '.rg.dups_marked.metrics.txt', rg_out)
    cmd = 'picard MarkDuplicates INPUT=' + rg_out + ' OUTPUT=' + mark_dup_out + ' METRICS_FILE=' + metrics_out
    subprocess.check_call(cmd, shell=True)



    # sam to bam and index
    cmd = 'samtools index ' + mark_dup_out
    subprocess.check_call(cmd, shell=True)

    os.remove(sam_out)
    os.remove(clean_out) 
    os.remove(sort_out)
    os.remove(rg_out)
