#!/bin/bash

#Created by Mihir Metkar on 01/29/16

#This script processes a .bam file after removal of XS field and only mapped alingments.
#Also provides some preliminary analysis about the bam file.

echo -e "Input gff subsetted itx file (itx__gff)"

read itx_file

echo -e "\nProcessing $itx_file ..."

bfile=`basename $itx_file`

fname="${bfile%.*}"

#total number of lines
echo -e "\nTotal number of lines in itx file..."
wc -l $itx_file | awk '{print $1}'

#distribution of type of entry- S, I or D
echo -e "\nNumber and type of entries..."
grep -v "^@" $itx_file | cut -f1 | sort | uniq -c | awk '{print $2 "\t" $1}'

#Calculate number of junctions in each read.
echo -e "\nCalculating number of Direct interactions for each read..."
echo -e "itx_type\t#_junctions\tCount"
awk '($1 == "D"){print $2}' $itx_file | uniq -c | awk '{print $2 "\t" $1}' | cut -f2 | sort -k1,1n | uniq -c | \
                awk '{print "D" "\t"$2 "\t" $1}'

#Calculate total type of interactions- Inter Vs Intra
echo -e "\nCalculating number of types of interactions..."
echo -e "itx_type\tCount"
awk '($1 == "D")' $itx_file | cut -f3,13 | awk '{if ($1 == $2) print "Intra-RNA"; else print "Inter-RNA"}' \
                | sort | uniq -c | awk '{print $2 "\t" $1}'

#Inter v Intra RNA ligations analysis
echo -e "\n\n"
read -t 10 -p "Do you want to do a Inter-RNA and Intra-RNA count analysis? (y/n) [default is 'y'] " answer
[ -z "$answer" ] && answer="y"  # if 'yes' have to be default choice
if [ "$answer" = "y" ]; then
  echo -e "\nCalculating number of Inter-RNA and Intra-RNA ligations...\n"
  echo -e "Type\tSpecies1\tSpecies2\tCount"

  awk '($1 == "D")' $itx_file | cut -f3,13 | awk '{if ($1 == $2) print "Intra-RNA""\t"$1"\t"$2; else print "Inter-RNA""\t"$1"\t"$2}' \
                  | sort -k1,1 -k2,2 -k3,3 | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'

else
  echo -e "\nSkipping counting Inter-RNA and Intra-RNA ligations\n"
fi

#Calculate number of singletons matching to each sequence in reference
echo -e "\n\n"
read -t 10 -p "Do you want to know how many singletons mapped to each reference sequence? (y/n) [default is 'y']" answer
[ -z "$answer" ] && answer="y"  # if 'yes' have to be default choice
if [ "$answer" = "y" ]; then
  echo -e "\nCalculating number of singletons mapping to reference...\n"
  echo -e "Species\tCount(singletons)"

  awk '($1 == "S"){print $3}' $itx_file | sort | uniq -c | awk '{print $2 "\t" $1}'

else
  echo -e "\nSkipping counting reference wise singleton counts\n"
fi


echo -e "\nDone\n"
