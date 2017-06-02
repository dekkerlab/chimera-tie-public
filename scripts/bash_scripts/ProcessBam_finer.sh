#!/bin/bash

#Created by Mihir Metkar on 01/29/16

#This script processes a .bam file after removal of XS field and only mapped alingments.
#Also provides some preliminary analysis about the bam file.

echo -e "Input bam file (*_noXS.bam)"

read bam_file

echo -e "\nProcessing $bam_file ..."

bfile=`basename $bam_file`

fname="${bfile%.*}"

#total number of reads (noXS file)
echo -e "\nTotal number of unqiue reads in bam file..."
samtools view $bam_file | wc -l

#Flag stats
echo -e "\nFLAG stats for bam file...\n"
echo -e "FLAG\tCount"
samtools view $bam_file | cut -f2 | sort | uniq -c | awk '{print $2 "\t" $1}'

#Get distribution of fragments numbers
echo -e "\nDistribution of fragment numbers from bam file...\n"
echo -e "#_fragments\tCount"
samtools view $bam_file | cut -f1 | uniq -c | awk '{print $1}' | sort -k1,1n | uniq -c | awk '{print $2 "\t" $1}'

#count number of reads with multiple fragments
echo -e "\nNumbers reads in bam file with multiple fragments i.e ligations..."
samtools view $bam_file | cut -f1 |  sort -k1,1n | uniq -cd | awk '{print $2}' | wc -l

#Get match length distribution
echo -e "\nPrinting match length distribution to file (length, count)...\n"
samtools view $bam_file | cut -f16 | cut -d ":" -f3 | sort -k1,1n | uniq -c | awk '{print $2 "\t" $1}' > ${fname}_matchLength.tsv

#Get distribution of number of reference matches. i.e how many reads matched to each reference element
#But asking a question here in case there are too many elements in reference
#if whiptail --yesno "Is this a good question" 2 6 ;then
read -t 10 -p "Do you want to know the number of fragments matching each reference sequence? (y/n) [default is 'y']" answer
[ -z "$answer" ] && answer="y"  # if 'yes' have to be default choice

if [ "$answer" = "y" ]; then
  echo -e "\nCalculating the number of reads matching to each reference element...\n"
  echo -e "reference\tcount"
  samtools view $bam_file | cut -f3 | sort | uniq -c | awk '{print $2 "\t" $1}'
else
  echo -e "\nSkipping counting # of matches to each reference!"
fi

echo -e "\nDone\n"
