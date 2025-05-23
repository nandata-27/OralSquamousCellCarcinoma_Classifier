#!/bin/bash
#SBATCH --job-name=Star_Align
#SBATCH --mail-type=ALL
#SBATCH --partition=standard
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G              
#SBATCH --time=1-00:00:00       
#SBATCH -o "outFile.txt"
#SBATCH -e "outErrorFile.txt"


#download fastq files 
module load sra-toolkit/3.2.0
file="/global/home/sa105231/Final_project/dl_fastq/SRR_Acc_List.txt"
download_dir="/global/home/sa105231/Final_project/dl_fastq"

while read -r srr; do  # Read each SRR ID line-by-line
    echo "Downloading $srr..."

    # Download the SRA file
    prefetch --max-size 100G "$srr" -O "$download_dir"

    # Convert to FASTQ (no .sra extension needed)
    fasterq-dump "$download_dir/$srr" -O "$download_dir"

    echo "Finished processing $srr"
done < "$file"

#clean fastq files - remove adapters, low quality reads
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

#align fastq files to reference human genome file 
#define directories
fastq_dir="/global/home/sa105231/Final_project/trimmed_fastq"
gen_ref="/global/home/sa105231/Final_project/GenomeReferenceDirectory/Genome"

#load module 
module load star/2.7.11b

#loop through all clean fastq files
for file in "${fastq_dir}"/*_1.trim.fastq; 
do

#Basename fuction removes the .read1.fastq ext and gets prefix part of fastq file without file path
    sample_name=$(basename "$file" _1.trim.fastq)
    if [ ! -f "${fastq_dir}/${sample_name}.sortedByCoord.out.bam" ]; then
        echo "Running aligningment on $file"
        #run star aligner
        STAR \
            --runThreadN 24 \
            --genomeDir ${gen_ref} \
            --readFilesIn ${sample_name}_1.trim.fastq ${sample_name}_2.trim.fastq \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.6 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMattributes NH HI nM MD \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${sample_name}
      else
        echo "Alignment performed for $file, skipping."
    fi
done

#script end
echo "finished aligning"

#index bam files 
#fastq directory
fastq_dir="/global/home/sa105231/Final_project/trimmed_fastq/processed_files/"

#load samtools
module load samtools

#generate BAM indexes

for file in "$fastq_dir"/*Aligned.sortedByCoord.out.bam; do
   samtools index ${file}
done

# generate counts
#fastq directory
fastq_dir="/global/home/sa105231/Final_project/trimmed_fastq/processed_files/"
gtf_file="/global/home/sa105231/Final_project/GenomeReferenceDirectory/Homo_sapiens.GRCh38.113.chr.gtf"
counts_dir="/global/home/sa105231/Final_project/trimmed_fastq/processed_files/Counts/"

#module load python and htseq
module load python/3.11
virtualenv htseq_env_py311
source htseq_env_py311/bin/activate
pip install htseq --no-index

#run htseq
for file in "$fastq_dir"/*.bam; do
sample_name=$(basename "$file" .bam)
echo "Processing file: $file"
    htseq-count -f bam -r pos -s reverse -m intersection-nonempty "$file" "$gtf_file" > "$counts_dir/${sample_name}_counts.txt"
done

#script end
echo "finished generating counts"



