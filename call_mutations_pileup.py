import sys
import os
import re
import fnmatch
import argparse
import subprocess

#loops through bam directory, runs mpileup and then calls the parsing script pileup_parser.py
#pilep_parser.py could be integrated as a function if so desired

#usage:
#module add python
#python call_mutations_pileup.py --bam_dir project_dir/bam/ --ref_fasta project_dir/reference_fasta.fa --out_dir project_dir/snv

#output directory is made by this script if it doesn't already exist

def parseArgs():

	parser = argparse.ArgumentParser()
	parser.add_argument('--bam_dir', required=True)
	parser.add_argument('--ref_fasta', required=True)
	parser.add_argument('--out_dir', default='.')

	args = vars(parser.parse_args())

	return args


if __name__ == '__main__':

    args = parseArgs()

    out_dir = args['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for f in os.listdir(args['bam_dir']):
        if fnmatch.fnmatch(f, '*.dups_marked.bam'): #if processing the bam steps differ, change to reflect
            splits = f.split('.')
            sample_name = splits[0]
            f = os.path.join(args['bam_dir'], f)

            sample_pileup = os.path.join(out_dir, sample_name + '.samtools.pileup')
            cmd = 'module add samtools/1.3.1; samtools mpileup -BQO -d 1000000 -q 20 -Q 20 -f' + args['ref_fasta']  + ' ' + f + ' >' + sample_pileup
            subprocess.check_call(cmd, shell=True)

            pileup_table = re.sub('.samtools.pileup', '.callstats.txt', sample_pileup)
            print('parsing ' + str(sample_pileup))
            cmd2 = 'python pileup_parser.py ' + sample_pileup + '> ' + pileup_table
            subprocess.check_call(cmd2, shell=True)
            #sys.exit()
                                    
          
