import sys
import re
import os
import fnmatch
import subprocess
import argparse

#generates flagstat and insertion size metrics for all bams in direcotyr

#usage: 
# module add python
# python qc_metrics_samples.py --ref project_dir/reference_fasta.fa --bam_dir project_dir/bam/ --out_dir project_dir/qc/
# qc folder is made if it doesn't exist already

def parseArgs():
  
    parser = argparse.ArgumentParser(description='Calculate summary target and sample coverage statistics.')
    #parser.add_argument('--gatk', default='/pkg/gatk/3.5/Linux/bin/GenomeAnalysisTK.jar')
    parser.add_argument('--ref', required=True)
    parser.add_argument('--bam_dir', required=True)
    #parser.add_argument('--min_base_quality', type=int, default=30);
    #parser.add_argument('--min_map_qual', type=int, default=15, help='Only count reads with at least this mapping quality (1)')
    parser.add_argument('--out_dir', required=True)

    args = vars(parser.parse_args())

    return args


def sampleName(bam):
    name = os.path.basename(bam)
    sample = name.split('.')[0]

    return sample



if __name__ == '__main__':

    args = parseArgs()

    out_dir = args['out_dir']
    #out_dir = os.path.join(args['out_dir'], 'qc')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    samples = []
    for file in os.listdir(args['bam_dir']):
        if fnmatch.fnmatch(file, '*.dups_marked.bam'):
            bam = os.path.join(args['bam_dir'], file)
            samples.append(sampleName(bam))
            cmd = 'module add R; module add java/jre/1.8.0_66; module add samtools/1.3.1; module add fastx; module add picard/2.6.0; python qc_metrics_sample.py' + ' --bam ' + bam + ' --ref ' + args['ref'] +  ' --out_dir ' + out_dir
            print(cmd)
            subprocess.check_call(cmd, shell=True)
	    #os.remove(os.path.join(out_dir, sampleName(bam) + '.sample_interval_summary'))
	    #os.remove(os.path.join(out_dir, sampleName(bam) + '.sample_interval_statistics'))
  

