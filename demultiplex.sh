#!/usr/bin/bash

# Preprocessing script for Gorilla project

##############################
#Demultiplexing Gorilla Data
##############################
cd /parfreylab/mann/gorillaYaws

#fastqc run
cp /parfreylab/shared/raw_data/miseq_data/18s/2019-March-ParfreyBilly_CustomLibraryG/ParfreyBilly_CustomLibraryG_L001_*gz .
gzip -d *gz
mkidr fastqc
fastqc ParfreyBilly_CustomLibraryG_L001_I1_001.fastq ParfreyBilly_CustomLibraryG_L001_I2_001.fastq ParfreyBilly_CustomLibraryG_L001_R1_001.fastq ParfreyBilly_CustomLibraryG_L001_R2_001.fastq -o fastqc
rm *fastq

#demultiplex
#this was done locally by Evan, if you want to run on the cluster need to install idemp and change path below
/projects/git/idemp/idemp -b demultiplexing_file_18s.revcomp.txt -I1 raw_data/ParfreyBilly_CustomLibraryG_L001_I1_001.fastq.gz -R1 raw_data/ParfreyBilly_CustomLibraryG_L001_R1_001.fastq.gz -R2 raw_data/ParfreyBilly_CustomLibraryG_L001_R2_001.fastq.gz -m 1 -o raw_data/demultiplexed/

#get some stats on the demultiplexed reads
cd demultiplex
ls *gz > file
ls *gz | while read line; do zcat $line | wc -l ; done >> raw
awk '{print $1/4}' raw > reads
paste file reads > readcount.txt
rm file raw reads
#get summary stats of raw reads
average=$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' readcount.txt)
median=$(awk '{print $2}' readcount.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else pri
nt a[x-1]; }')
max=$(awk '{print $2}' readcount.txt | sort -n | tail -n1)
min=$(awk '{print $2}' readcount.txt | sort -n | head -n1)
echo "Average:" $average >> readcount.stats.txt
echo "Median:" $median >> readcount.stats.txt
echo "Minimum read count:" $min >> readcount.stats.txt
echo "Maximum read count:" $max >> readcount.stats.txt

#copy demultiplexed data to entamoeba, ready for dada2 processing