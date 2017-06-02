#!/bin/bash

#Created by Mihir Metkar on 11/22/16.

#This script takes in itx2scatter_strand.py input and makes a Pymol distance file,
#to measure and draw distance lines in Pymol

echo -e "Input scatter plot file (Output of itx2scatter_strand.py)"
read spfile

echo -e "\nInput PDB molecule name"
read pdb_name

echo -e "\nCut-off for interaction count. [x >= cutoff](Use 0 if you want all interactions)"
read cutoff

echo -e "\nProcessing file $spfile ...\n"

fname="${spfile%.*}"

while IFS=$'\t' read -r col1 col2 col3 col4 col5
do
  if [ "$col5" -ge "$cutoff" ]
  then
  #dist test, 4ug0_5S and resi 2 and name P, 4ug0_5S and resi 119 and name P
  echo "dist ${col1}_${col2}_$cutoff, $pdb_name and resi $col3 and name P, $pdb_name and resi $col4 and name P" >> ${fname}_${cutoff}.pml
  fi
done < "$spfile"

echo -e "\nDone\n"
