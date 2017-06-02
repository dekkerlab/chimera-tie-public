#!/bin/bash

#Created by Mihir Metkar on 03/27/17

#This script processes matrices to create multiple scaling plots.
echo -e "Following script processes multiple unzipped matrices
(required format- *.matrix) to make scaing plot text file.
It will extract binsize from filename (format- *_10nt_*.matrix).\n"

#For each file in folder, process it
for matrix in *.matrix;
do
  file=$matrix
  infile=`basename $matrix`
  fname="${infile%.*}"
  echo -e "Input file:\n$infile"
  #Split file name by "_" and pull out binsize by matching pattern.
  binsize=$(echo "$fname" | gawk -F_ '{ for (i = 1; i <= NF; ++i) if ($i ~ /nt/) print $i } ' \
            | sed 's/[^0-9]*//g')
  echo -e "binsize is $binsize"
  perl ~/scripts/Scaling2_batch.pl $matrix $binsize

  echo -e "\n";
done
