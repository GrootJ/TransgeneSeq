#!/usr/bin/env pythport os
import sys
import subprocess
import os
import re
import os.path
import argparse
import glob
#import tempfileimport glob
import shutil
#import tarfile

#usage: via qc_metrics_samples.py

def parseArgs():
  
    parser = argparse.ArgumentParser(description='Calculate summary target and sample coverage statistics.')
    #parser.add_argument('--gatk', default='/pkg/gatk/3.5/Linux/bin/GenomeAnalysisTK.jar')
    parser.add_argument('--ref', required=True)
    parser.add_argument('--bam', required=True)
    #parser.add_argument('--min_base_quality', type=int, default=30)
    #parser.add_argument('--min_map_qual', type=int, default=15, help='Only count reads with at least this mapping quality (1)')
    parser.add_argument('--out_dir', default='.')

    args = vars(parser.parse_args())

    return args


def sampleName(bam):
    name = os.path.basename(bam)
    sample = name.split('.')[0]

    return sample


def generate_ref_index(ref):

    fai_file = ref + '.fai'
    if not os.path.exists(fai_file):
        cmd = 'samtools faidx ' + ref # generate .fai file
        subprocess.check_call(cmd, shell=True)

    dict_file = re.sub('.fa', '.dict', ref)
    if not os.path.exists(dict_file):
        cmd = 'picard CreateSequenceDictionary R=' + ref + ' O=' + dict_file 
        subprocess.check_call(cmd, shell=True) 


def flagstat(bam, ref, out_dir):
    outFile = os.path.join(out_dir, sampleName(bam) + ".flagstat")
    cmd = '/camhpc/pkg/gatk/3.5/centos6/bin/GenomeAnalysisTK -T FlagStat --filter_bases_not_stored -R ' + ref + ' -I ' + bam + ' -o ' + outFile
    subprocess.check_call(cmd, shell=True)

def insert_size(bam, out_dir):
    outFile = os.path.join(out_dir, sampleName(bam) + '.insert_size_metrics.txt')
    hist = os.path.join(out_dir, sampleName(bam) + '.histogram.pdf')
    cmd = 'picard CollectInsertSizeMetrics H='+ hist + ' I=' + bam + ' O=' + outFile
    subprocess.check_call(cmd, shell=True)
#def depth_of_coverage(bam, gatk, ref, min_base_quality, min_map_qual, out_dir):
    #param_prefix = os.path.join(out_dir, sampleName(bam))
    #cmd = 'java -Xmx15g -jar ' + gatk + ' --downsampling_type NONE' + ' --validation_strictness SILENT' + ' -R ' + ref + ' -T DepthOfCoverage -o ' + param_prefix + ' -I ' + bam + ' --interval_merging OVERLAPPING_ONLY --minBaseQuality ' + str(min_base_quality) + ' --minMappingQuality ' + str(min_map_qual) + ' --omitLocusTable --omitPerSampleStats --countType COUNT_FRAGMENTS'
    #subprocess.check_call(cmd, shell=True)
    #if os.path.exists(param_prefix): os.rename(param_prefix, param_prefix + '.sample_base_summary')


if __name__ == '__main__':
  
    args = parseArgs()
    generate_ref_index(args['ref'])
    flagstat(args['bam'], args['ref'], args['out_dir'])
    insert_size(args['bam'], args['out_dir'])
    #depth_of_coverage(args['bam'], args['gatk'], args['ref'], args['min_base_quality'], args['min_map_qual'], args['out_dir'])
    

