#!/bin/bash 
#SBATCH --job-name=removeAdapter
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24qg15@queensu.ca
#SBATCH --partition=standard
#SBATCH -c 2
#SBATCH --mem=128G             
#SBATCH --time=1-00:00:00       
#SBATCH -o "output.txt"
#SBATCH -e "error.txt"


#file paths
adapter="/global/home/sa105231/Final_project/adapters.fa"
fastq_dir="/global/home/sa105231/Final_project/"

#load module 
module load bbmap/39.06

#run bbduk on all fastq files
for file in "${fastq_dir}"/*_1.fastq; 
do
    #Basename fuction removes the .read1.fastq ext and gets prefix part of fastq file without file path
    sample_name=$(basename "$file" _1.fastq)
    
    #run bbduk 
    bbduk.sh \
    in1=${sample_name}_1.fastq in2=${sample_name}_2.fastq \
    out1=${sample_name}_1.trim.fastq out2=${sample_name}_2.trim.fastq \
    ref=$adapter \
    k=21 \
    tpe=t \
    ktrim=r \
    useshortkmers=f \
    mink=7 \
    qtrim=rl \
    trimq=15 \
    minlength=30 \
    trimpolygleft=5 trimpolygright=5 \
    lhist=${sample_name}.lhist.txt stats=${sample_name}.stats.txt

done

echo "Finished adapter trimming and contaminant filtering"