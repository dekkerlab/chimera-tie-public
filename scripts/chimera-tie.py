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
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    
    args=parser.parse_args()

    fastq_1=args.fastq_1
    fastq_2=args.fastq_2
    bowtie2_idx_prefix=args.bowtie2_idx_prefix
    verbose=args.verbose
    
    side1_file_name=os.path.basename(fastq_1)
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
    sam_file=side1_file_name+'___'+genome_name+'.sam'
    output = open(sam_file, "w")

    bowtie_cmd = bowtie_path+ ' --local -p 8 -x '+bowtie2_idx_prefix+' -U '+fastq_1+' -S '+sam_file
    verboseprint(bowtie_cmd)
    bowtie_args = shlex.split(bowtie_cmd)
    verboseprint(bowtie_args)
    bowtie_proc = subprocess.Popen(bowtie_args,
                        stdout=output,
                        stderr=subprocess.PIPE,)

    err = bowtie_proc.communicate()
    bowtie_proc.wait()
    
    sam_fh=input_wrapper(sam_file)
    
    for i,line in enumerate(sam_fh):
        line=line.rstrip("\n")
        if line.startswith("#") or line.startswith("@"):
            continue
        
        x=line.split("\t")
        
        #unmapped or secondary
        #if( (int(x[1]) & 0x4) or (int(x[1]) & 0x100) ):
        #   continue
            
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
        cigar_tup=zip(values[1::2],(values[0::2]))
        for i in cigar_tup:
            cigar_dict[i[0]]+=int(i[1])
        
        read_length=cigar_dict['M']+cigar_dict['I']+cigar_dict['S']+cigar_dict['=']+cigar_dict['X'];
        alignment_length=cigar_dict['M']+cigar_dict['=']+cigar_dict['X']+cigar_dict['D']+cigar_dict['N'];
        
        print(qname,flag,rname,pos,cigar,tlen,read_length,alignment_length,seq,len(seq))
        
    sam_fh.close
        

    
    

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
    
def output_wrapper(outfile):
    
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
    
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

if __name__=="__main__":
    main()

   