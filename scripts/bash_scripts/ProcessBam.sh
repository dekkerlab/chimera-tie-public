#!/bin/bash

#Created by Mihir Metkar on 11/22/16

#This script processes a chimera-tie.py output .bam file. Also provides some preliminary analysis about the bam file.

echo -e "Input bam file (output of chimera-tie.py)"

read bam_file

echo -e "\nProcessing $bam_file ..."

bfile=`basename $bam_file`

fname="${bfile%.*}"

# Get header from bam file
samtools view -H $bam_file > $fname.header
echo -e "\nWrote header file"

# Remove multi-mappers and unmapped reads
echo -e "\nPulling out uniquely mapped reads...\n"
samtools view $bam_file | grep -v "XS:i" | awk '($2 != 4)' | cat $fname.header - | samtools view -bS - > ${fname}_mapped_noXS.bam

echo -e "\nWrote bam file with uniquely mapped reads"

# Remove reads with more than 3 mismatches
bamtools filter -tag 'XM':'<=3' -in ${fname}_mapped_noXS.bam -out ${fname}_mapped_noXS_XM3.bam

echo -e "\nWrote bam file with maximum 3 mismatches"

# Remove reads with more than 3 gaps
bamtools filter -tag 'XG':'<=3' -in ${fname}_mapped_noXS_XM3.bam -out ${fname}_mapped_noXS_XM3XG3.bam

echo -e "\nWrote bam file with maximum 3 gaps"

# Count number of lines
read -t 10 -p "\nDo you want to count the number of lines? (y/n) [default is 'n']\n" answer
[ -z "$answer" ] && answer="n"  # if 'yes' have to be default choice

if [ "$answer" = "y" ]; then
  echo -e "\nNumber of lines in original bam file:"
  samtools view $bam_file | wc -l
else
  echo -e "\nSkipping counting number of lines!"
fi

#count sam flags
echo -e "\nsam flags and count:"
samtools view $bam_file | cut -f2 | sort | uniq -c

# Count lines in new file
read -t 10 -p "\nDo you want to count the number of lines in new file? (y/n) [default is 'n']\n" answer
[ -z "$answer" ] && answer="n"  # if 'no' have to be default choice

if [ "$answer" = "y" ]; then
  echo -e "\nNumber of lines in new bam file:"
  samtools view ${fname}_mapped_noXS_XM3XG3.bam | wc -l

else
  echo -e "\nSkipping counting number of lines in new bam file!"
fi

#count sam flags in new file
echo -e "\nsam flags and count:"
samtools view ${fname}_mapped_noXS.bam | cut -f2 | sort | uniq -c

#Count number of unique read IDs
echo -e "\nNumber of unique read IDs in the new bam file:"
samtools view ${fname}_mapped_noXS.bam | cut -f1 | uniq | wc -l

echo -e "\nDone\n"
