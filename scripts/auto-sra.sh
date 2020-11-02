#!/bin/bash

#take user input for specific BioProject number
echo Please enter the NCBI BioProject number: 
read bp_num

#keeping things nice and tidy
echo Please enter the directory to which all files will be transferred:
read down_dir
mkdir ~/projects/bioinfo/seq-dl/$down_dir
dir=~/projects/bioinfo/seq-dl/$down_dir

#fetching and downloading metadata for given bioproject number
echo Obtaining metadata on project $bp_num...
esearch -db sra -query PRJNA257197 | efetch -format runinfo > $dir/runinfo.csv

#update
echo Output metadata on project $bp_num as $dir/runinfo.csv.

#extract all SRR numbers in associated with bioproject number
echo Extracting SRR numbers from metadata...
cat $dir/runinfo.csv | cut -f 1 -d ',' | grep SRR > $dir/runids.txt
echo Extraction of SRR numbers complete. Output save as $dir/runids.txt

#download fastq
echo Obtaining fastq files...
cat $dir/runids.txt | parallel "echo Downloading {} ..."
cat $dir/runids.txt | parallel fasterq-dump --split-files {} -O $dir/fastq -e 12 -t ~/dev/shm -p -x
echo Download complete. All files located at $dir
