from __future__ import print_function
from __future__ import division

import numpy as np
import sys
import argparse
import time
import re
import gzip
import os
import math

#samtools view file.bam | sort -V -k1,1 -k13,13 | python RipitClash_chimeratie.py -ref exonco-ords.txt > youroutputname.itx

def main():
	# Get input from user. Look up function argparse (argument parse)
	parser=argparse.ArgumentParser(description='RipitClash pipline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-i',help='input file for RipitClash',dest='infile',default=None,type=str)
	#parser.add_argument('-ref',help='reference file for chrome',dest='ref_file',default=None,required=True,type=str)
	#parser.add_argument('-b',help='binsize in bp',dest='binsize',default=100,type=int)

	args=parser.parse_args()

	infile=args.infile

	#if infile is zipped, it unzips the file
	if infile==None:
		infh=sys.stdin
	else:
		if infile.endswith('.gz'):
			infh=gzip.open(infile,'r')
		else:
			infh=open(infile,'r')

	cis=0
	trans=0
	direct=0
	indirect=0
	gap=0
	chrs=dict()
	type=None
	junction=None

	#initialize last line
	#Better way to do define an empty list?
	lastLineList=[None,None,None,None,None,None,None,None,None,None,None,None,None,None]
	lastLine=None
	# Split all columns in the input file by tab and put it in a [list] file
	for lineNumber,line in enumerate(infh): #enumerate gives a number to each item.
		currentLine=line.rstrip("\n") #?
		currentLineList=currentLine.split("\t") #split by tab and make a list
		#print (lastLineList, currentLineList)

		#define items in currentLineList
		current_readID=currentLineList[0] #read ID, field 1
		current_chrID=currentLineList[2] #chr ID, field 3
		currentLineList[3]=int(currentLineList[3]) #start of fragment, position in reference, convert it to integer
		current_start=currentLineList[3]
		#where do we ask it to store previous line in lastLine?
		last_readID=lastLineList[0]
		last_chrID=lastLineList[2]
		#lastLineList[3]=int(lastLineList[3]) #Cant do this as first value is "None"
		last_start=lastLineList[3]

		#what is this doing? Count number of interactions on each chromosome
		if current_chrID in chrs:
			chrs[current_chrID] += 1
		else:
			chrs[current_chrID] = 1

		#To identify fragments within a read, pull out fragments only if readIDs are the same
		#classify interaction as cis or trans, right now w.r.t each chromosome, so somewhat irrelevant
		if last_readID == current_readID: #for a particular read
			if current_chrID == last_chrID: #if chr ID is the same
				cis = cis + 1
				type="cis"
			else:
				trans = trans + 1
				type="trans"

			#Define START and END in read for each fragment
			#Since >12 fields are only present in mapped reads, cant define it earlier?
			last_xx=lastLineList[12]
			last_xy=lastLineList[13]
			last_xx=int(last_xx.split(":")[2])
			last_xy=int(last_xy.split(":")[2])

			current_xx=currentLineList[12]
			current_xy=currentLineList[13]
			current_xx=int(current_xx.split(":")[2]) #split by : and take 2nd filed
			current_xy=int(current_xy.split(":")[2])

			if last_xy==current_xx:
				#print (last_readID, last_xx, last_xy, current_readID, current_xx, current_xy,sep="\t")
				junction='direct'
				direct+=1

			else:
				#print (last_readID, last_xx, last_xy, current_readID, current_xx, current_xy,sep="\t")
				gap=current_xx-last_xy
				junction='indirect'
				indirect+=1

			#check if file is sorted. If read IDs dont match or XQ numbers are not sorted, it will give an error.
			if ((lastLineList[0] != currentLineList[0])or (current_xx <= last_xx)):
				print ("Use -V to sort")
				sys.exit("ERROR: file out of order! cur_xx="+str(current_xx)+" < ,last_xx="+str(last_xx))

		#else:
			#print (last_readID, current_readID)
			#print ("starting new read!")

			#print the smaller value first, may not be a good idea.
			#Because then we miss the direction of ligation between different ligations in the same read

			print(type,last_readID,last_chrID,last_start,current_readID,current_chrID,current_start,junction,"gap="+str(gap),sep="\t")

			#if we decide to assign ligation point to positions- 3' of of acceptor and 5' end of donor
			ligation_pt_last=last_start+last_xy
			ligation_pt_current=current_start

			print(type,last_readID,last_chrID,ligation_pt_last,current_readID,current_chrID,ligation_pt_current,junction,"gap="+str(gap),sep="\t")
			print ("\n")
			
			gap=0
		lastLine = currentLine
		lastLineList = currentLineList
	sys.stderr.write ('Direct ligations='+str(direct)+'\n')
	sys.stderr.write ('Indirect ligations='+str(indirect)+'\n')
	sys.stderr.write ('Intra chromosome interactions='+str(cis)+'\n')
	sys.stderr.write ('Inter chromosome interactions='+str(trans)+'\n')


	infh.close()

if __name__=="__main__":
	main()
