#!/usr/bin/env python

# encoding: utf-8

"""
chimeraTie suite created by Bryan R. Lajoie on 09/16/2015
This script created by Mihir Metkar on 07/19/2016
"""

from __future__ import print_function
from __future__ import division

import sys
import argparse
import subprocess
import shlex
import logging
import time
import re
import os
import math
import uuid
import socket
from collections import defaultdict
from collections import Counter
from datetime import datetime

verboseprint=lambda *a, **k: None
__version__ = "1.0"

def main():

    parser=argparse.ArgumentParser(description='Modify Emsembl gtf, output of DEXseq_anotation.py, to chimeraTie standards',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-g', '--gene_annotation', dest='gene_annotation', type=str, required=True, help='path to gene annoation GFF file')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    gene_annotation=args.gene_annotation
    verbose=args.verbose

    verbose=args.verbose

    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)

    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None

    verboseprint("\n",end="")

    # load GFF file
    verboseprint("processing GFF file ... ")
    gtf_open=open(gene_annotation,"r")

    annotation_file_name=get_file_name(gene_annotation)

    outfile_transcripts=annotation_file_name+"_transcripts.gff"
    outfile_exons=annotation_file_name+"_exons.gff"

    transcripts_out_fh=open(outfile_transcripts,'w+')

    print ("#bin","chrom","processed","entry_type","txStart","txEnd","strand","name","exonic_part_number","name2",sep="\t",file=transcripts_out_fh)

    exons_out_fh=open(outfile_exons,'w+')

    print ("#bin","chrom","processed","entry_type","txStart","txEnd","strand","name","exonic_part_number","transcript_id","name2",sep="\t",file=exons_out_fh)


    for index,line in enumerate(gtf_open):
        x=line.rstrip("\n")
        x=line.split("\t")

        """
        1	dexseq_prepare_annotation.py	aggregate_gene	11869	14412	.	+	.	gene_id "ENSG00000223972"
        1	dexseq_prepare_annotation.py	exonic_part	11869	11871	.	+	.	transcripts "ENST00000456328"; exonic_part_number "001"; gene_id "ENSG00000223972"
        """

        chrom=x[0]
        processed=x[1]
        entry_type=x[2]
        txStart=x[3]
        txEnd=x[4]
        dot=x[5]=x[7]
        strand=x[6]

        if processed!="dexseq_prepare_annotation.py":
            print ("Wrong input! Convert gtf to gff using dexseq_prepare_annotation.py")
            break

        if entry_type=="aggregate_gene":
            name=x[8].split("\"")[1]
            exonic_part_number=0
            name2=name+"_"+str(exonic_part_number)

            print (index,chrom,processed,entry_type,txStart,txEnd,strand,name,exonic_part_number,name2,sep="\t",file=transcripts_out_fh)

        elif entry_type=="exonic_part":
            #transcripts "ENST00000456328"; exonic_part_number "001"; gene_id "ENSG00000223972"

            transcript_id=x[8].split("\"")[1]
            exonic_part_number=x[8].split("\"")[3]
            name=x[8].split("\"")[5]
            name2=name+"_"+exonic_part_number


            print (index,chrom,processed,entry_type,txStart,txEnd,strand,name,exonic_part_number,transcript_id,name2,sep="\t",file=exons_out_fh)


    transcripts_out_fh.close()
    exons_out_fh.close()
    gtf_open.close()

def get_file_name(file):

    file_name=file.split("/")[-1]

    short_name = file_name
    file_name = re.sub(r"\.matrix\.gz$", "", file_name)
    file_name = re.sub(r"\.matrix$", "", file_name)
    file_name = re.sub(r"\.gz$", "", file_name)

    # if non-matrix file - remove extension
    if short_name == file_name:
        file_name=remove_file_extension(file_name)

    return file_name

def remove_file_extension(file):

    tmp = file.split(".")
    del tmp[-1]
    file_name = '.'.join(tmp)

    return file_name

if __name__=="__main__":
    main()
