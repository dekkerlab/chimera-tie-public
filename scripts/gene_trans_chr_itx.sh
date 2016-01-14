#!/bin/bash
#For a particular chromosome subsetted itx file, this scripts extracts trans interactions for a particular gene.
gene=${1} #gene name
itx_file=${2} #location of the sorted & subsetted itx file
output_file_prefix=${3}
echo ""
echo "Extracting trans interactions for ${gene} from ${itx_file}."
echo ""
echo "*Ensure itx file was subsetted before using this script.*"
echo ""

grep "${gene}" ${itx_file} | awk "((\$19==\"${gene}\") && (\$23!=\"${gene}\"))||((\$19!=\"${gene}\") && (\$23==\"${gene}\")){print \$0}" > ${gene}_${output_file_prefix}_trans.itx
echo " lines per file"
echo "----------------"
wc -l ${gene}_${output_file_prefix}_trans.itx
awk "(\$23==\"${gene}\"){print \$12}" ${gene}_${output_file_prefix}_trans.itx > loc1_${gene}
wc -l loc1_${gene}
awk "(\$19==\"${gene}\"){print \$4}" ${gene}_${output_file_prefix}_trans.itx > loc2_${gene}
wc -l loc2_${gene}
cat loc2_${gene} loc1_${gene} > loc12_${gene}
wc -l loc12_${gene}
sort -n loc12_${gene} | uniq -c > loc_${gene}_count
wc -l loc_${gene}_count
awk "{print \$1 \"\t\" \$2}" loc_${gene}_count > loc_${gene}_count_${output_file_prefix}.tsv
wc -l loc_${gene}_count_${output_file_prefix}.tsv
rm loc1_${gene} && rm loc2_${gene} && rm loc12_${gene} && rm loc_${gene}_count
echo ""

echo "done"
