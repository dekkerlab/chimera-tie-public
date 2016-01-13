#Classifying intra-chromosome, inter-RNA interactions.
"""Usage: (output sorted by starting position of location 1 and 2)
python classify_itx.py --trans -i test.pairwise-score.txt | sort -k2,2n -k5,5n > test_output.tsv"""
from __future__ import print_function
from __future__ import division
from collections import defaultdict
from collections import Counter
from operator import itemgetter
from collections import OrderedDict
import sys
import argparse
import re
import gzip
import os
import math

parser=argparse.ArgumentParser(description='classify chr wide interactions',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i',help='input file, text format',dest='infile',default=None,type=str)
parser.add_argument('-o',help='output file name, sam format',dest='outfile',default="output.txt",type=str)
parser.add_argument('-t',help='threshold for itx score',dest='threshold',default="10",type=int)
parser.add_argument('--trans',help='only get trans itx',dest='trans',action='store_true') #store_true= false!?

args=parser.parse_args()
trans=args.trans
#print (trans)
threshold=args.threshold

outfile=args.outfile
#ofh=open(outfile,'w')

infile=args.infile
if infile==None:
	infh=sys.stdin
else:
	if infile.endswith('.gz'):
		infh=gzip.open(infile,'r')
	else:
		infh=open(infile,'r')

def main():
	Total_itx=0
	gene1_itx_counter=[]
	for line in infh:
		if line.startswith("yHeader"):
			continue
		column = line.split('\t')
		gene1 = column[0]
		gene2 = column[1]
		if column[2]=="NA": #Is this right? What is NA?
			score=0
		score = int(float(column[2]))
		Total_itx+=1
		if score < threshold:
			continue
		gene1_itx_counter.append(gene1)

		gene1_name=gene1.split('|')[0]
		gene1_loc=gene1.split('|')[2]
		gene1_start=int(re.split('[:-]',gene1_loc)[1])
		gene1_end=int(re.split('[:-]',gene1_loc)[2])
		gene2_name=gene2.split('|')[0]
		gene2_loc=gene2.split('|')[2]
		gene2_start=int(re.split('[:-]',gene2_loc)[1])
		gene2_end=int(re.split('[:-]',gene2_loc)[2])
		#print (gene1, gene1_loc,gene1_start, gene1_end,gene2, gene2_loc,gene2_start,gene1_end, sep="\t")
		if trans==True and gene1==gene2:
			continue
		else:
			dist=abs(gene1_end-gene2_start)
			print (gene1_name,gene1_start, gene1_end, gene2_name,gene2_start,gene2_end, score, dist, sep="\t")


	print ("Total interactions on chr:",Total_itx, sep="\t")
	print ("Total number of genes:",len(gene1_itx_counter), sep="\t")
	num_itx(gene1_itx_counter)
	#print (num_itx_col1, sep="\t")

	infh.close()
	#ofh.close()

def num_itx(gene_list):
	#number of times this gene is present in the itx file above the threshold itx
	#print (gene_list)
	count=Counter(gene_list)
	sorted_count = OrderedDict(sorted(count.items(), key=itemgetter(1)))
	for gene,cnt in sorted_count.iteritems():
		print (gene, cnt, sep="\t")
			#return (gene, cnt)

if __name__=="__main__":
	main()
