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

bin_dir = sys.path[0] + "/"

def main():

    parser=argparse.ArgumentParser(description='A split-read alignment process for mapping multi-chimeric reads',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f1', '--fastq1', dest='fastq_1', type=str, required=True, help='input side 1 fastq file')
    parser.add_argument('-f2', '--fastq2', dest='fastq_2', type=str, required=True, help='input side 2 fastq file')
    parser.add_argument('-g', '--genome', dest='bowtie2_idx_prefix', type=str, required=True, help='path to bowtie2 idx file')
    parser.add_argument('-l', '--localmode', dest='bowtie2_local', type=str, default='sensitive-local', required=True, help='local alignment mode')
    parser.add_argument('-b', '--bopts', dest='bowtie2_options', nargs='+', type=str, default=[], required=False, help='optional bowtie2 alignment options')
    parser.add_argument('--minseq', dest='min_seq_len', type=int, default=10, help='minimum sequence length to attempt alignment')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    
    args=parser.parse_args()

    fastq_1=args.fastq_1
    fastq_2=args.fastq_2
    min_seq_len=args.min_seq_len
    bowtie2_local=args.bowtie2_local
    bowtie2_options=args.bowtie2_options
    bowtie2_idx_prefix=args.bowtie2_idx_prefix
    verbose=args.verbose
    
    genome_name=os.path.basename(bowtie2_idx_prefix)
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    
    verboseprint("\n",end="")

    # check that bowtie2 and samtools are installed
    check_bowtie()
    check_samtools()
    
    # check that the bowtie2 idx exists
    check_bowtie_index(bowtie2_idx_prefix)
    
    # perform alignment (side1)  
    
    # perform initial cleanup if necessary
    fastq_1_name=os.path.basename(fastq_1)
    sam_file=fastq_1_name+'___'+genome_name+'.sam'
    try:
        os.remove(sam_file)
    except OSError:
        pass
    
    fastq_file=fastq_1
    i=0
    while os.stat(fastq_file).st_size > 0:
        print("iteration",i)
        fastq_file=iterate_chimera_tie(i,fastq_file,sam_file,bowtie2_idx_prefix,genome_name,bowtie_path,min_seq_len,bowtie2_local,bowtie2_options)
        i+=1
    

def iterate_chimera_tie(iter,fastq_file,sam_file,bowtie2_idx_prefix,genome_name,bowtie_path,min_seq_len,bowtie2_local,bowtie2_options):

    sam_fh = open(sam_file, "a")
    
    print(iter,fastq_file)
    
    # open same file
    bowtie_sam_file=genome_name+'.bowtie.sam'
    bowtie_sam_fh = open(bowtie_sam_file, "w")

    bowtie2_options_str=' '.join(bowtie2_options)

    bowtie_cmd = bowtie_path+ ' '+'--'+bowtie2_local+' '+bowtie2_options_str+' -x '+bowtie2_idx_prefix+' -U '+fastq_file+' -S '+bowtie_sam_file
    print(bowtie_cmd)

    bowtie_args = shlex.split(bowtie_cmd)
    bowtie_proc = subprocess.Popen(bowtie_args,
                        stdout=bowtie_sam_fh,)
    bowtie_proc.wait()
    
    bowtie_sam_fh.close()    
        
    # open fastq file
    fastq_file=genome_name+'__'+str(iter)+'.fastq'
    fastq_fh = open(fastq_file, "w")
    
    # read sam file
    bowtie_sam_fh=input_wrapper(bowtie_sam_file)
    
    for i,line in enumerate(bowtie_sam_fh):
        line=line.rstrip("\n")
        if line.startswith("#") or line.startswith("@"):
            if iter == 0:
                print(line,file=sam_fh)
            continue
        
        x=line.split("\t")
        
        qname=x[0]
        flag=x[1]
        rname=x[2]
        pos=x[3]
        mapq=x[4]
        cigar=x[5]
        tlen=x[8]
        seq=x[9]
        qual=x[10]
        
        cigar_dict = Counter({'M':0, 'I':0, 'D':0, 'N':0, 'S':0, 'H':0, 'P':0, 'X':0, '=':0})
        
        pattern = re.compile('([MIDNSHPX=])')
        values = pattern.split(cigar)[:-1]
        cigar_tup=zip(values[0::2],values[1::2])

        # reverse cigar if reverse strand
        if(int(x[1]) & 0x10):
            cigar_tup=cigar_tup[::-1]
            flat = [x for sublist in cigar_tup for x in sublist]
            cigar=''.join(flat)
            qual=qual[::-1]
        
        for i in cigar_tup:
            cigar_dict[i[1]]+=int(i[0])
        
        # check later to ensure this is correct
        read_length=cigar_dict['M']+cigar_dict['I']+cigar_dict['S']+cigar_dict['=']+cigar_dict['X']
        match_length=cigar_dict['M']+cigar_dict['=']+cigar_dict['X']+cigar_dict['N']+cigar_dict['I']
        span_length=cigar_dict['M']+cigar_dict['=']+cigar_dict['X']+cigar_dict['N']+cigar_dict['D']
        
        offset=None
        if(len(qname.split(":::")) == 2):
            offset=qname.split(":::")[-1]
        qname=qname.split(":::")[0]
        offset_start=offset_end=0
        if(offset != None):
            offset_start,offset_end=offset.split("-")
  
        xx=int(offset_start)
        
        if(len(cigar_tup) > 1):
            left_cigar=cigar_tup[0]
            right_cigar=cigar_tup[-1]    
            
            if(left_cigar[1] == 'S'):
                left_seq=seq[0:int(left_cigar[0])]
                left_qual=qual[0:int(left_cigar[0])]
                left_start=int(offset_start)
                left_end=int(left_cigar[0])+int(offset_start)
                xx += int(left_cigar[0])

                if(len(left_seq) > min_seq_len):
                    print("@"+qname+":::"+str(left_start)+"-"+str(left_end),left_seq,"+",left_qual,sep="\n",file=fastq_fh)

            if(right_cigar[1] == 'S'):
                right_seq=seq[len(seq)-int(right_cigar[0]):len(seq)]
                right_qual=qual[len(qual)-int(right_cigar[0]):len(qual)]
                right_start=len(seq)-int(right_cigar[0])+int(offset_start)
                right_end=len(seq)+int(offset_start)

                if(len(right_seq) > min_seq_len):
                    print("@"+qname+":::"+str(right_start)+"-"+str(right_end),right_seq,"+",right_qual,sep="\n",file=fastq_fh)
        
        xy=xx+match_length

        # capture every line to aggregrate sam
        if cigar != "*":
            tmp=line.split("\t")
            tmp.insert(12,"ZR:i:"+str(read_length))
            tmp.insert(12,"ZS:i:"+str(span_length))
            tmp.insert(12,"ZM:i:"+str(match_length))
            tmp.insert(12,"XY:i:"+str(xy))
            tmp.insert(12,"XX:i:"+str(xx))
            line='\t'.join(tmp)

        print(line,file=sam_fh)

    bowtie_sam_fh.close()
    fastq_fh.close()
    
    sam_fh.close()
    
    return(fastq_file)
    
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
    
def sam_fh_wrapper(outfile):
    
    if outfile.endswith('.gz'):
        fh=gzip.open(outfile,'wb')
    else:
        fh=open(outfile,'w')
    
    suppress_comments=0
    
    # disable comment(s)if (UCSC format file)
    if outfile.endswith('.bed'):
        suppress_comments = 1
    if outfile.endswith('.bed.gz'):
        suppress_comments = 1
    if outfile.endswith('.bedGraph'):
        suppress_comments = 1
    if outfile.endswith('.bedGraph.gz'):
        suppress_comments = 1
    if outfile.endswith('.wig'):
        suppress_comments = 1
    if outfile.endswith('.wig.gz'):
        suppress_comments = 1

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

   
   