#!/bin/bash

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
