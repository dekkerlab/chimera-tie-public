#!/bin/bash

#Created by Mihir Metkar on 03/07/17

#This script processes minus ligase and ligase matrices to make a dual heatmap.

echo -e "\nInput directory for minus ligase matrices (without last /)"
read minus_dir
echo -e "\nInput directory for plus ligase matrices (without last /)"
read plus_dir
echo -e "\nInput txt file with gene names, usually the D interaction count file"
read count_file
echo -e "\nprefix for output files"
read outfile_name

cat $count_file | while read transcript count
do
  echo -e "\nregion_name : $transcript"
  #for format ngene_id | transcript_id | transcript_name, use following cmd
  tname=$(echo $transcript | awk -F\| '{print $3}')

  #For the name as is
  #tname=$transcript
  
  echo -e "\n$tname"
  minus_matrix=$(ls $minus_dir/ | grep ${tname} | grep 'matrix.gz$')
  plus_matrix=$(ls $plus_dir/ | grep ${tname} | grep 'matrix.gz$')

  echo -e "\n$minus_matrix\n$plus_matrix"

  infile=`basename $plus_matrix`
  fname="${infile%%.*}"

  perl ~/bin/github/cworld-dekker/scripts/perl/heatmap.pl \
  -i $minus_dir/$minus_matrix -i $plus_dir/$plus_matrix --em --dl -o ${outfile_name}_${tname}

  echo -e "\nWrote matrix!"

done
