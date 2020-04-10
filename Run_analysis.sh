#!/usr/bin/env bash

#run example: Run_analysis.sh projectname projectPath fastqPath referencePath

#project name ; example:"ASYN_univector_GS_pYZZ002_Asyn-gs"
project_name=$0
#project path; example:"/home/dlin2/Projects/CHO/TransgeneSeq/Train2/"
project_path=$1
#fastq_path; example:"/camhpc/ngs/TechDev/Transgene_Seq/test_new_pipeline/training_2019/FASTQ/fastq/"
fastq_path=$2
#reference_path; example:"/camhpc/ngs/TechDev/Transgene_Seq/Transgene_References/transgene_vectors/"
reference_path=$3

mkdir ${project_path}/FASTQ
mkdir ${project_path}/prealign_qc
mkdir ${project_path}/postalign_qc
mkdir ${project_path}/R_analysis
mkdir ${project_path}/bam
mkdir ${project_path}/stdout
mkdir ${project_path}/snv
mkdir ${project_path}/reference

cp ${fastq_path}/*.fastq.gz ../FASTQ/
ln -s ${reference_path}/${project_name}.fa ../reference/${project_name}.fa

#qc
fastqc --noextract ../FASTQ/*.fastq.gz  -o ../prealign_qc/
##multiqc ../prealign_qc/ -o ../prealign_qc/


module load python/2.7.11

#alignment
python align_samples_qsub.py --fastq_dir ${project_path}/FASTQ/ --ref_fasta ${project_path}/reference/${project_name}.fa --out_dir ${project_path}/bam/  --stdout_dir ${project_path}/stdout/

#Picard QC
python qc_metrics_samples.py --ref ${project_path}/reference/${project_name}.fa  --bam_dir ${project_path}/bam/ --out_dir ${project_path}/postalign_qc/

#mutation call
python call_mutations_pileup.py --bam_dir ${project_path}/bam/  --ref_fasta ${project_path}/reference/${project_name}.fa  --out_dir ${project_path}/snv/

# Run Rscript
module load R

