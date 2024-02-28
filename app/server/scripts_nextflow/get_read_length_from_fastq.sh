#!/usr/bin/env bash


###### Function to extract read length from fastq files 
echo "--------------------------------------------------"
echo "#################### READ LENGTHS ################"
echo "--------------------------------------------------"

fastq_file="$1" # $input_dir/*/*/*/*.fastq.gz
sample="$2" # ERRxyz
output_dir="$3"

#Make necessary environments and files
mkdir -p $output_dir
if [ ! -f $output_dir/"$sample"_read_lengths_pass.txt ]
then
    echo "Length" > $output_dir/"$sample"_read_lengths_pass.txt
fi

#Filter read length information from fastq files
if [[ "$fastq_file" == *".gz" ]]
then
    echo "HERE" 
    zcat $fastq_file | awk 'NR%4==2' | awk '{ print length }' >> $output_dir/"$sample"_read_lengths_pass.txt
else
    echo "HERE2"
    cat $fastq_file | awk 'NR%4==2' | awk '{ print length }' >> $output_dir/"$sample"_read_lengths_pass.txt
fi

