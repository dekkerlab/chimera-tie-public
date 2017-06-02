#!/bin/bash

#Created by Mihir Metkar on 03/07/17

#This script processes subsetted itx file to make multiple matrices.

echo -e "\nInput gff subsetted itx file (itx__gff)"
read itx_file
echo -e "\nInput bam header file"
read header
echo -e "\nInput cutoff for number of D interactions for each gene. (< #D > cutoff will be selected. Eg. 0,400)"
read cutoff

echo -e "\nProcessing $itx_file with $cutoff as cutoff...\n"

infile=`basename $itx_file`
fname="${infile%.*}"

count_file=${fname}_transcripts_D_${cutoff}.txt

awk '($1 == "D")' $itx_file | cut -f3 | sort | uniq -c | sort -k1,1nr | \
awk -v cut=$cutoff '($1 > cut)' | awk '{print $2 "\t" $1}' > $count_file

echo -e "\nWrote a file with Direct interaction count for each gene. (More than cutoff)...\n"

echo -e "\n\n"
read -t 10 -p "Are you processing minus ligase (y/n) [default is 'n']" answer
[ -z "$answer" ] && answer="n"  # if 'yes' have to be default choice
if [ "$answer" = "y" ]; then
  echo -e "\n\tProvide a Direct interaction count file of Plus ligase...\n"
  read count_file

else
  echo -e "\n\tUsing the recently created Direct interaction count file.\n"

fi

cat $count_file | while read transcript count
do
  echo "region_name : $transcript, count : $count"

  region=\'$transcript\' #method 2 to get region in the format geneID\|trans_id\|name --> region="${transcript//|/\\|}"
  tname=$(grep "${transcript}" $header | awk '{print $2}' | awk -F\| '{print $3}')
  length=$(grep "${transcript}" $header | awk '{print $3}' | awk -F: '{print $2}')
  bsize=$(awk -v len=$length 'BEGIN{print int(len/100)}')

  echo -e "gene_name : $tname, length : $length, bsize : $bsize\n"

  writeMatrix="python ~/bin/github/chimera-tie/scripts/itx2matrix_hakan_sym.py -i $itx_file -r $region -D --bsize $bsize -o ${fname}_${tname}_${bsize}nt_D.matrix.gz"

  eval $writeMatrix
done

echo -e "\nWrote all matrices with approximately 100 bins\n"
