#!/usr/bin/env python

# encoding: utf-8

"""
chimera-tie.py
Created by Bryan R. Lajoie on 09/16/2015
Some base functions taken from tophat
https://github.com/infphilo/tophat
-
Created by Cole Trapnell on 2008-12-25.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
Updated and maintained by Daehwan Kim and Geo Pertea since Jul 2010.
-
"""

from __future__ import print_function
from __future__ import division

import numpy as np
import scipy as sp
import sys
import argparse
import subprocess
import shlex
import logging
import itertools
import time
import gzip
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
debug = None

bin_dir = sys.path[0] + "/"

#############################################################################

def get_arguments():

    parser=argparse.ArgumentParser(description='Construct interaction matrix from chimeraTie ITX file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--itx', dest='itx_file', type=str, required=True, help='input chimeraTie ITX file')
    parser.add_argument('-r', '--regions', dest='regions', type=str, required=False,
                         help='Regions to be picked from the interaction file')
    parser.add_argument('--bsize', dest='bin_size', type=int, default=1000)
    parser.add_argument('-S', dest='singleton', action='store_true', help='include singletons in matrix')
    parser.add_argument('-I', dest='indirect', action='store_true', help='include indirect itx in the matrix')
    parser.add_argument('-D', dest='direct', action='store_true', help='include direct itx in the matrix')
    parser.add_argument('-g', dest='genome', type=str, required=False, default = 'hg19',  help='Genome Assembly')
    parser.add_argument('--debug', dest='debug', action='store_true', help='debug mode')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    return parser.parse_args()

#############################################################################

def get_region_list(region_argument):
    if region_argument is None:
        return None
    regions = region_argument.split(",")
    if len(regions) == 0:
        return None
    return regions

#############################################################################

def get_headers(itx_file, bin_size, genome, regions):

    genes=dict()
    n_bins=0
    header_rows=[]
    header_cols=[]

    GENE_NAME_INDEX = 0
    CHR_INDEX = 1
    START_INDEX = 3
    LENGTH_INDEX = 4

    itx_name=get_file_name(itx_file)

    itx_fh=input_wrapper(itx_file)
    for line_num,line in enumerate(itx_fh):
        verboseprint(line)
        x=line.rstrip("\n").split("\t")

        if line.startswith("#"):
            continue

        if line.startswith("@"):
            x=line.lstrip("@").rstrip("\n").split("\t")

            if x[GENE_NAME_INDEX] not in regions:
                continue

            n_gene_bins=int(math.ceil(int(x[LENGTH_INDEX])/bin_size))
            gene_start=int(x[START_INDEX])
            gene_length=int(x[LENGTH_INDEX])
            gene_end=gene_start+gene_length

            genes[x[GENE_NAME_INDEX]]=n_bins

            for i in xrange(n_gene_bins):
                bin_start=int(x[START_INDEX])+(i*bin_size)
                bin_end=bin_start+bin_size
                if bin_end > gene_end:
                    bin_end = gene_end

                header=str(x[GENE_NAME_INDEX])+'__'+str(i)+'|' +\
                         str(genome) +'|'+\
                         str(x[CHR_INDEX])+':'+str(bin_start)+'-'+str(bin_end)
                header_rows.append(header)
                header_cols.append(header)

            n_bins+=n_gene_bins

            continue

        break

    itx_fh.close()

    if n_bins > 50000:
            sys.exit('danger - matrix would be too large!, use itx2subset.py first')

    return { 'genes': genes,
             'header_rows': header_rows,
             'header_cols': header_cols,
             'n_bins': n_bins }


#############################################################################

def arrange_logging(verbose):
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)

def is_circular( frag1_start, frag2_start,
                 frag1_xx, frag2_xx,
                 frag1_strand, frag2_strand ):

    circular = False

    if (frag1_xx > frag2_xx and \
       frag1_strand == '+' and frag2_strand == '+') or \
       (frag1_xx < frag2_xx and \
       frag1_strand == '-' and frag2_strand == '-') or \
       (frag1_xx < frag2_xx and \
       frag1_strand == '+' and frag2_strand == '-') or \
       (frag1_xx > frag2_xx and \
       frag1_strand == '-' and frag2_strand == '+'):
       circular = True

    return circular

#############################################################################

def get_matrix(itx_file, n_bins, bin_size, genes,
               singleton, indirect, direct, regions):
    itx_fh=input_wrapper(itx_file)

    matrix=np.zeros([n_bins,n_bins])
    matrix_sum=0

    REGION_1_NAME_INDEX = 22
    REGION_2_NAME_INDEX = 27

    FRAG1_CHROM_INDEX = 2
    FRAG1_START_INDEX = 3
    FRAG1_STRAND_INDEX = 5
    FRAG1_MATCHLENGTH_INDEX = 9
    FRAG1_GENENAME_INDEX = REGION_1_NAME_INDEX
    FRAG1_GENE_CHROM_INDEX = 23
    FRAG1_GENE_START_INDEX = 25
    FRAG1_GENE_END_INDEX = 26
    FRAG1_XX_INDEX = 6

    FRAG2_CHROM_INDEX = 12
    FRAG2_START_INDEX = 13
    FRAG2_STRAND_INDEX = 15
    FRAG2_MATCHLENGTH_INDEX = 19
    FRAG2_GENENAME_INDEX = REGION_2_NAME_INDEX
    FRAG2_GENE_CHROM_INDEX = 26
    FRAG2_GENE_START_INDEX = 30
    FRAG2_GENE_END_INDEX = 31
    FRAG2_XX_INDEX = 16


    for line_num, line in enumerate(itx_fh):
        x=line.rstrip("\n").split("\t")

        if line.startswith("#"):
            continue

        if line.startswith("@"):
            continue

        if len(x) < 32:
            continue

        if x[REGION_1_NAME_INDEX] not in regions or\
           x[REGION_2_NAME_INDEX] not in regions:
           continue

        interaction_type=x[0]

        # only keep specified itx by type
        if not singleton and interaction_type == 'S':
            continue
        if not indirect and interaction_type == 'I':
            continue
        if not direct and interaction_type == 'D':
            continue

        frag1_chrom=x[FRAG1_CHROM_INDEX]
        frag1_start=int(x[FRAG1_START_INDEX])
        frag1_matchlength=int((x[FRAG1_MATCHLENGTH_INDEX].split(":")[-1]))
        frag1_end=frag1_start+frag1_matchlength
        frag1_gene_name=x[FRAG1_GENENAME_INDEX]
        frag1_gene_chrom=x[FRAG1_GENE_CHROM_INDEX]
        frag1_gene_start=int(x[FRAG1_GENE_START_INDEX])
        frag1_gene_end=int(x[FRAG1_GENE_END_INDEX])

        frag2_chrom=x[FRAG2_CHROM_INDEX]
        frag2_start=int(x[FRAG2_START_INDEX])
        frag2_matchlength=int((x[FRAG2_MATCHLENGTH_INDEX].split(":")[-1]))
        frag2_end=frag2_start+frag2_matchlength
        frag2_gene_name=x[FRAG2_GENENAME_INDEX]
        frag2_gene_chrom=x[FRAG2_GENE_CHROM_INDEX]
        frag2_gene_start=int(x[FRAG2_GENE_START_INDEX])
        frag2_gene_end=int(x[FRAG2_GENE_END_INDEX])

        #print(line_num,interaction_type,frag1_chrom,frag1_start,frag1_matchlength,frag1_end, sep="\t")
        #print(frag1_gene_name,frag1_gene_chrom,frag1_gene_start,frag1_gene_end,frag2_gene_name,frag2_gene_chrom,frag2_gene_start,frag2_gene_end, sep="\t")

        frag1_txn_coord_start_offset=frag1_start-frag1_gene_start
        frag1_txn_coord_end_offset=frag1_end-frag1_gene_start
        frag2_txn_coord_start_offset=frag2_start-frag2_gene_start
        frag2_txn_coord_end_offset=frag2_end-frag2_gene_start

        #print(frag1_trn_coord_start,frag1_trn_coord_end,frag2_trn_coord_start,frag2_trn_coord_end, sep="\t")

        frag1_txn_bin_start=int(math.floor(frag1_txn_coord_start_offset/bin_size))+genes[frag1_gene_name]
        frag1_txn_bin_end=int(math.floor(frag1_txn_coord_end_offset/bin_size))+genes[frag1_gene_name]+1
        frag2_txn_bin_start=int(math.floor(frag2_txn_coord_start_offset/bin_size))+genes[frag2_gene_name]
        frag2_txn_bin_end=int(math.floor(frag2_txn_coord_end_offset/bin_size))+genes[frag2_gene_name]+1

        #print (frag1_txn_coord_start_offset,frag1_txn_coord_end_offset,frag2_txn_coord_start_offset,frag2_txn_coord_end_offset, frag1_txn_bin_start,frag1_txn_bin_end,frag2_txn_bin_start,frag2_txn_bin_end,sep="\t")
        frag1_strand = x[FRAG1_STRAND_INDEX]
        frag2_strand = x[FRAG2_STRAND_INDEX]

        matrix_sum += 1

        if interaction_type == 'S':
            frag1_interacting_nucleotide = (frag1_txn_bin_start + \
                                            frag1_txn_bin_end) / 2
            frag2_interacting_nucleotide = frag1_interacting_nucleotide
            matrix[frag1_interacting_nucleotide][frag2_interacting_nucleotide] += 1
            continue

        # We have a non-singleton case
        frag1_xx = int(x[FRAG1_XX_INDEX].split(":")[-1] )
        frag2_xx = int(x[FRAG2_XX_INDEX].split(":")[-1] )
        circular = is_circular( frag1_start, frag2_start,
                                frag1_xx, frag2_xx,
                                frag1_strand, frag2_strand )


        if circular == True:
            if frag1_strand == '+':
                frag1_interacting_nucleotide = frag1_txn_bin_end
            else:
                frag1_interacting_nucleotide = frag1_txn_bin_start

            if frag2_strand == '-':
                frag2_interacting_nucleotide = frag2_txn_bin_end
            else:
                frag2_interacting_nucleotide = frag2_txn_bin_start

        matrix[frag1_interacting_nucleotide][frag2_interacting_nucleotide] += 1
        matrix[frag2_interacting_nucleotide][frag1_interacting_nucleotide] += 1
        #matrix[frag1_txn_bin_start:frag1_txn_bin_end,:][:,frag2_txn_bin_start:frag2_txn_bin_end] += 1
        #matrix[frag2_txn_bin_start:frag2_txn_bin_end,:][:,frag1_txn_bin_start:frag1_txn_bin_end] += 1



    itx_fh.close()
    verboseprint("\twrote",matrix_sum,"itx")
    return matrix

#############################################################################

def main():

    args      = get_arguments()
    itx_file  = args.itx_file
    bin_size  = args.bin_size
    regions   = get_region_list(args.regions)
    global debug
    debug     = args.debug
    verbose   = args.verbose

    arrange_logging(verbose)

    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None

    headers = get_headers(itx_file, bin_size, args.genome, regions)
    n_bins = headers['n_bins']

    verboseprint(n_bins,"x",n_bins)

    itx_name=get_file_name(itx_file)
    matrix = get_matrix(itx_file,
                         n_bins = n_bins,
                         bin_size = bin_size,
                         genes = headers['genes'],
                         singleton=args.singleton,
                         indirect = args.indirect,
                         direct = args.direct,
                         regions = regions)

    print("itx name is ", itx_name, "\n\n")
    matrixFile=itx_name+".matrix.gz"
    writeMatrix(headers['header_rows'], headers['header_cols'],
                matrix, matrixFile)

############################################################

def writeMatrix(header_rows,header_cols,matrix,matrixFile,precision=4):
    """
    write a np matrix with row/col headers - my5C file format - txt formatted gzipped file
    """

    nrows=len(header_rows)
    ncols=len(header_cols)

    # interaction matrix output
    out_fh=gzip.open(matrixFile,"wb")

    # write matrix col headers
    header=[str(i) for i in header_cols]
    print(str(nrows)+"x"+str(ncols)+"\t"+"\t".join(header),file=out_fh)

    format_func=("{:0."+str(precision)+"f}").format

    k=0

    for i in xrange(nrows):
        print(header_rows[i]+"\t"+"\t".join(map(format_func,matrix[i,:])),file=out_fh)

    out_fh.close

############################################################

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

############################################################

def remove_file_extension(file):

    tmp = file.split(".")
    del tmp[-1]

    file_name = '.'.join(tmp)

    return file_name

############################################################

def which(program):
    def is_executable(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        progpath = os.path.join(bin_dir, program)
        if is_executable(progpath):
           return progpath
        for path in os.environ["PATH"].split(os.pathsep):
           progpath = os.path.join(path, program)
           if is_executable(progpath):
              return progpath
    return None

############################################################

def die(msg=None):
  if msg is not None:
    print(msg)
  sys.exit(1)

############################################################

def prog_path(program):
    progpath=which(program)
    if progpath == None:
        die("Error locating program: "+program)
    return progpath

############################################################

def input_wrapper(infile):
    if infile.endswith('.gz'):
        fh=gzip.open(infile,'r')
    else:
        fh=open(infile,'r')

    return fh

############################################################

def output_wrapper(outfile,append=False,suppress_comments=False):

    if outfile.endswith('.gz'):
        if append:
            fh=gzip.open(outfile,'a')
        else:
            fh=gzip.open(outfile,'w')
    else:
        if append:
            fh=open(outfile,'a')
        else:
            fh=open(outfile,'w')

    # disable comment(s)if (UCSC format file)
    if outfile.endswith('.bed'):
        suppress_comments = True
    if outfile.endswith('.bed.gz'):
        suppress_comments = True
    if outfile.endswith('.bedGraph'):
        suppress_comments = True
    if outfile.endswith('.bedGraph.gz'):
        suppress_comments = True
    if outfile.endswith('.wig'):
        suppress_comments = True
    if outfile.endswith('.wig.gz'):
        suppress_comments = True
    if outfile.endswith('.sam'):
        suppress_comments = True
    if outfile.endswith('.sam.gz'):
        suppress_comments = True
    if outfile.endswith('.bam'):
        suppress_comments = True
    if outfile.endswith('.bam.gz'):
        suppress_comments = True

    if not suppress_comments:
        print("## ",os.path.basename(__file__),sep="",file=fh)
        print("## ",sep="",file=fh)
        print("## Dekker Lab",sep="",file=fh)
        print("## Contact: Bryan R. Lajoie",sep="",file=fh)
        print("## https://github.com/blajoie",sep="",file=fh)
        print("## ",sep="",file=fh)
        print("## Version:\t",__version__,sep="",file=fh)
        print("## Date:\t",get_date(),sep="",file=fh)
        print("## Host:\t",get_compute_resource(),sep="",file=fh)

    return(fh)


#############################################################

def get_date():
    time=datetime.now()
    date=time.strftime('%I:%M:%S %p, %m/%d/%Y')

    return date

#############################################################

def get_compute_resource():
    return(socket.gethostname())

#############################################################

if __name__=="__main__":
    main()
