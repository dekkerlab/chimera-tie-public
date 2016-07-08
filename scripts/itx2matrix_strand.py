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

def main():

    parser=argparse.ArgumentParser(description='Construct interaction matrix from chimeraTie ITX file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--itx', dest='itx_file', type=str, required=True, help='input chimeraTie ITX file')
    parser.add_argument('--bsize', dest='bin_size', type=int, default=1000)
    parser.add_argument('-S', dest='singleton', action='store_true', help='include singletons in matrix')
    parser.add_argument('-I', dest='indirect', action='store_true', help='include indirect itx in the matrix')
    parser.add_argument('-D', dest='direct', action='store_true', help='include direct itx in the matrix')
    parser.add_argument('--debug', dest='debug', action='store_true', help='debug mode')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    itx_file=args.itx_file
    bin_size=args.bin_size
    singleton=args.singleton
    indirect=args.indirect
    direct=args.direct
    global debug
    debug=args.debug
    verbose=args.verbose

    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)

    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None

    genes=dict()
    n_bins=0
    header_rows=[]
    header_cols=[]

    itx_name=get_file_name(itx_file)

    itx_fh=input_wrapper(itx_file)
    for line_num,line in enumerate(itx_fh):
        x=line.rstrip("\n").split("\t")

        if line.startswith("#"):
            continue

        if line.startswith("@"):
            x=line.lstrip("@").rstrip("\n").split("\t")

            n_gene_bins=int(math.ceil(int(x[4])/bin_size))
            gene_start=int(x[3])
            gene_length=int(x[4])
            gene_end=gene_start+gene_length

            genes[x[0]]=n_bins

            for i in xrange(n_gene_bins):
                bin_start=int(x[3])+(i*bin_size) #xrange generates intergers, 0 to n. Here, 0*bin_size, 1*bin_size etc
                bin_end=bin_start+bin_size
                if bin_end > gene_end:
                    bin_end = gene_end

                header=str(x[0])+'__'+str(i)+'|hg19|'+str(x[1])+':'+str(bin_start)+'-'+str(bin_end)
                header_rows.append(header)
                header_cols.append(header)

            n_bins+=n_gene_bins

            continue

        break

    itx_fh.close()

    if n_bins > 50000:
            sys.exit('danger - matrix would be too large!, use itx2subset.py first')

    verboseprint(n_bins,"x",n_bins)

    matrix=np.zeros([n_bins,n_bins])
    matrix_sum=0

    itx_name=get_file_name(itx_file)

    itx_fh=input_wrapper(itx_file)
    for line_num,line in enumerate(itx_fh):
        x=line.rstrip("\n").split("\t")

        if line.startswith("#"):
            continue

        if line.startswith("@"):
            continue

        interaction_type=x[0]

        # only keep specified itx by type
        if not singleton and interaction_type == 'S':
            continue
        if not indirect and interaction_type == 'I':
            continue
        if not direct and interaction_type == 'D':
            continue

        frag1_chrom=x[2]
        frag1_start=int(x[3])
        frag1_matchlength=int((x[9].split(":")[-1]))
        frag1_end=frag1_start+frag1_matchlength
        frag1_gene_name=x[22]
        frag1_gene_chrom=x[23]
        frag1_gene_start=int(x[25])
        frag1_gene_end=int(x[26])

        frag2_chrom=x[12]
        frag2_start=int(x[13])
        frag2_matchlength=int((x[19].split(":")[-1]))
        frag2_end=frag2_start+frag2_matchlength
        frag2_gene_name=x[27]
        frag2_gene_chrom=x[28]
        frag2_gene_start=int(x[30])
        frag2_gene_end=int(x[31])

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

        matrix[frag1_txn_bin_start:frag1_txn_bin_end,:][:,frag2_txn_bin_start:frag2_txn_bin_end] += 1
        matrix[frag2_txn_bin_start:frag2_txn_bin_end,:][:,frag1_txn_bin_start:frag1_txn_bin_end] += 1

        matrix_sum += 1

    verboseprint("\twrote",matrix_sum,"itx")

    matrixFile=itx_name+".matrix.gz"
    writeMatrix(header_rows,header_cols,matrix,matrixFile)


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

def load_gff(gene_annotation):

    if not os.path.isfile(gene_annotation):
        die("GFF file is missing"+gene_annotation)

    genes=dict()
    ignored_genes=set()

    header2index = dict()

    gff_fh=input_wrapper(gene_annotation)

    for i,line in enumerate(gff_fh):
        line=line.rstrip("\n")
        x=line.split("\t")

        if line.startswith("#"):
            # header
            #bin name    chrom   strand  txStart txEnd   cdsStart    cdsEnd  exonCount   exonStarts  exonEnds    score   name2   cdsStartStat    cdsEndStat  exonFrames
            for index,header in enumerate(x):
                header2index[header]=index
            continue

        # process GFF entiries
        name=x[header2index["name"]]
        name2=x[header2index["name2"]]
        chrom=x[header2index["chrom"]]
        # remove any haplotype info
        chrom=chrom.split("_")[0]
        start=x[header2index["txStart"]]
        end=x[header2index["txEnd"]]

        if name2 in ignored_genes:
            continue
        elif name2 not in genes:
            genes[name2]=defaultdict(list)
            genes[name2]["chrom"]=chrom
            genes[name2]["start"]=start
            genes[name2]["end"]=end
        else:
            if chrom != genes[name2]["chrom"]:
                verboseprint("WARNING: ignoring "+name2+" (exists on multiiple chromosomes!)")
                ignored_genes.add(name2)
                del genes[name2]
                continue

            if start < genes[name2]["start"]:
                genes[name2]["start"]=start
            if end > genes[name2]["end"]:
               genes[name2]["end"]=end

    n_genes=len(genes)
    n_ignored_genes=len(ignored_genes)
    verboseprint("\tignored",n_ignored_genes," genes in GFF")
    verboseprint("\tfound",n_genes,"genes in GFF")
    verboseprint("")

    return(genes)

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

def count_lines(file):
    fh = input_wrapper(file)

    count = 0
    for _ in fh:
        count += 1
    return count

def revcomp(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    dna = dna[::-1]
    result = [complement.get(base,"N") for base in dna]
    return ''.join(result)

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

def die(msg=None):
  if msg is not None:
    print(msg)
  sys.exit(1)

def prog_path(program):
    progpath=which(program)
    if progpath == None:
        die("Error locating program: "+program)
    return progpath

def check_bowtie_index(idx_prefix, add="(genome)"):

    idxext="bt2"
    bowtie_ver="2 "

    idx_fwd_1 = idx_prefix + ".1."+idxext
    idx_fwd_2 = idx_prefix + ".2."+idxext
    idx_rev_1 = idx_prefix + ".rev.1."+idxext
    idx_rev_2 = idx_prefix + ".rev.2."+idxext

    if os.path.exists(idx_fwd_1) and \
       os.path.exists(idx_fwd_2) and \
       os.path.exists(idx_rev_1) and \
       os.path.exists(idx_rev_2):
        return

    bwtidxerr="Error: Could not find Bowtie "+bowtie_ver+"index files (" + idx_prefix + ".*."+idxext+")"

    bwtidx_env = os.environ.get("BOWTIE2_INDEXES")

    if bwtidx_env == None:
        die(bwtidxerr)
    if os.path.exists(bwtidx_env+idx_fwd_1) and \
       os.path.exists(bwtidx_env+idx_fwd_2) and \
       os.path.exists(bwtidx_env+idx_rev_1) and \
       os.path.exists(bwtidx_env+idx_rev_2):
        if os.path.exists(bwtidx_env + idx_prefix + ".1.ebwt") and os.path.exists(bwtidx_env + idx_prefix + ".1.bt2"):
            print >> sys.stderr, bwtbotherr
        return
    else:
        die(bwtidxerr)

def check_bowtie():
    bowtie_bin = "bowtie2"

    global bowtie_path
    bowtie_version = None
    bowtie_path=which(bowtie_bin)
    if bowtie_path:
      bowtie_version = get_bowtie_version()
    if bowtie_version == None:
           die("Error: Bowtie not found on this system.")

def check_samtools():
    samtools_bin = "samtools"

    global samtools_path
    samtools_version = None
    samtools_path=which(samtools_bin)
    if samtools_path == None:
           die("Error: samtools not found on this system.")

def get_bowtie_version():
    try:
        # Launch Bowtie to capture its version info
        proc = subprocess.Popen([bowtie_path, "--version"],
                          stdout=subprocess.PIPE)

        stdout_value = proc.communicate()[0]

        bowtie_version = None
        if not stdout_value: stdout_value=''
        bowtie_out = stdout_value.splitlines()[0]
        version_str=" version "
        ver_str_idx = bowtie_out.find(version_str)
        if ver_str_idx != -1:
            version_val = bowtie_out[(ver_str_idx + len(version_str)):]
            bvers=re.findall(r'\d+', version_val)
            bowtie_version = [int(x) for x in bvers]
        while len(bowtie_version)<4:
            bowtie_version.append(0)
        return bowtie_version
    except OSError, o:
       errmsg=fail_str+str(o)+"\n"
       if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
           errmsg+="Error: bowtie not found on this system"
       die(errmsg)

def input_wrapper(infile):
    if infile.endswith('.gz'):
        fh=gzip.open(infile,'r')
    else:
        fh=open(infile,'r')

    return fh

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

def get_date():
    time=datetime.now()
    date=time.strftime('%I:%M:%S %p, %m/%d/%Y')

    return date

def get_compute_resource():
    return(socket.gethostname())

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def getSmallUniqueString():
    tmp_uniq=str(uuid.uuid4())
    tmp_uniq=tmp_uniq.split('-')[-1]
    return(tmp_uniq)

def de_dupe_list(input):
    """de-dupe a list, preserving order.
    """

    sam_fh = []
    for x in input:
        if x not in sam_fh:
            sam_fh.append(x)
    return sam_fh

if __name__=="__main__":
    main()
