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

    parser=argparse.ArgumentParser(description='A split-read alignment process for mapping multi-chimeric reads',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f1', '--fastq1', dest='fastq_1', type=str, required=True, help='input side 1 fastq file')
    parser.add_argument('-f2', '--fastq2', dest='fastq_2', type=str, required=True, help='input side 2 fastq file')
    parser.add_argument('-x', '--bowtie_idx', dest='bowtie2_idx_prefix', type=str, required=True, help='path to bowtie2 idx file')
    parser.add_argument('-g', '--gene_annotaiton', dest='gene_annotation', type=str, default='path to gene annoation GFF file')
    parser.add_argument('--dd', '--distannce_definition', dest='distance_definition', type=int, default=0)
    parser.add_argument('--bopts', '--bopts', dest='bowtie2_options', type=str, default='--local -D 20 -R 3 -N 0 -L 16 -i S,1,0.50', required=False, help='optional bowtie2 alignment options')
    parser.add_argument('--debug', dest='debug', action='store_true', help='debug mode - keep itetation fastq/sam')
    parser.add_argument('--minseq', dest='min_seq_len', type=int, default=10, help='minimum sequence length to attempt alignment')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    
    args=parser.parse_args()

    fastq_1=args.fastq_1
    fastq_2=args.fastq_2
    min_seq_len=args.min_seq_len
    bowtie2_idx_prefix=args.bowtie2_idx_prefix
    gene_annotation=args.gene_annotation
    distance_definition=args.distance_definition
    bowtie2_options=args.bowtie2_options
    global debug
    debug=args.debug
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
        
    # perform initial cleanup if necessary
    fastq_1_name=get_file_name(fastq_1)
    header_file=fastq_1_name+'___'+genome_name+'.sam.header'
    unsorted_sam_file=fastq_1_name+'___'+genome_name+'.unsorted.sam'
    try:
        os.remove(unsorted_sam_file)
    except OSError:
        pass
    
    fastq_file=fastq_1
    
    verboseprint("processing GFF file ... ")
    load_gff(gene_annotation)

    verboseprint("running iterative-chimera-split ... \n")
    # run iterative/split mapping
    i=0
    while count_lines(fastq_file) > 0:
        print("iteration #",i,sep="")
        fastq_file=iterate_chimera_tie(i,fastq_file,unsorted_sam_file,header_file,bowtie2_idx_prefix,genome_name,bowtie_path,min_seq_len,bowtie2_options)
        i+=1

    if not debug:
        os.remove(fastq_file)

    verboseprint("sorting sam ... ")
    sorted_sam_file=fastq_1_name+'___'+genome_name+'.sorted.sam'
    sam_in_fh=input_wrapper(unsorted_sam_file)
    sam_out_fh=output_wrapper(sorted_sam_file)
    sort_cmd = "sort -V -k1,1 -k13,13"
    sort_args = shlex.split(sort_cmd)
    sort_proc = subprocess.Popen(sort_args,
                    stdin=sam_in_fh,
                    stdout=sam_out_fh,)
    sort_proc.wait()
    sam_in_fh.close()
    sam_out_fh.close()
    verboseprint("\tdone")

    # remove unsorted file
    if not debug:
        os.remove(unsorted_sam_file)

    verboseprint("")

    filenames = [header_file, sorted_sam_file]
    sam_file=fastq_1_name+'___'+genome_name+'.sam'
    with output_wrapper(sam_file) as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    if not debug:
        os.remove(sorted_sam_file)
        os.remove(header_file)

    verboseprint("writing bam ... ")
    samtools_cmd = "samtools view -bS "+sam_file
    bam_file=fastq_1_name+'___'+genome_name+'.bam'
    bam_fh=output_wrapper(bam_file)
    samtools_args = shlex.split(samtools_cmd)
    samtools_proc = subprocess.Popen(samtools_args,
                    stdout=bam_fh,)
    samtools_proc.wait()
    bam_fh.close()

    os.remove(sam_file)

    verboseprint("bam2itx ... ")
    bam2itx_cmd = "samtools view "+bam_file
    bam2itx_args = shlex.split(bam2itx_cmd)
    bam2itx_proc = subprocess.Popen(bam2itx_args,
                    stdout=subprocess.PIPE)
    
    nrow=500
    ncol=350
    matrix=np.zeros([nrow,ncol])
    matrix.fill(np.nan)
    header_rows=[]
    header_cols=range(ncol)
    matrix_index=0
    
    previous_read_id=None
    buffer=[]
    with bam2itx_proc.stdout:
        for line_num,line in enumerate(iter(bam2itx_proc.stdout.readline, b'')):
            x=line.split("\t")
            
            #unmapped or secondary
            if( (int(x[1]) & 0x4) or (int(x[1]) & 0x100) ):
                continue
            
            current_read_id=x[0]
            
            if((current_read_id != previous_read_id) and (previous_read_id != None) and (len(buffer) != 0)):
            
                for i1,b1 in enumerate(buffer):
                    tmp_b1=b1.split("\t")
                    
                    xx_1=int(tmp_b1[12].split(":")[-1])
                    xy_1=int(tmp_b1[13].split(":")[-1])
                    
                    for i2,b2 in enumerate(buffer):
                        if i1 >= i2:
                            continue
                            
                        tmp_b2=b2.split("\t")    
                        xx_2=int(tmp_b2[12].split(":")[-1])
                        xy_2=int(tmp_b2[13].split(":")[-1])
                        
                        i_type="I"
                        if abs(xy_1-xx_2) <= distance_definition:
                            i_type="D"
                    
                        #print(i_type,previous_read_id,tmp_b1[2],tmp_b1[3],tmp_b1[12],tmp_b1[13],tmp_b1[14],tmp_b1[15],tmp_b1[16],tmp_b1[17],tmp_b2[2],tmp_b2[3],tmp_b2[12],tmp_b2[13],tmp_b2[14],tmp_b2[15],tmp_b2[16],tmp_b2[17],sep="\t")
                
                if(matrix_index < nrow) and len(buffer) > 1:
                    start=int(buffer[0].split("\t")[12].split(":")[-1])
                    end=int(buffer[-1].split("\t")[13].split(":")[-1])
                    if(end < 350 and end > 300):
                        header_rows.append(previous_read_id)
                        for i,b in enumerate(buffer):
                            tmp_b=b.split("\t")
                            tmp_start=int(tmp_b[12].split(":")[-1])
                            tmp_end=int(tmp_b[13].split(":")[-1])
                            
                            score=0
                            if i % 2 == 0:
                                score=1
                            matrix[matrix_index,tmp_start:tmp_end]=score
                            
                        matrix_index += 1
                        
                buffer=[]
            
            # add current line to buffer 
            buffer.append(line)
            
            # set previous = current
            previous_read_id=current_read_id
    
        for i1,b1 in enumerate(buffer):
            tmp_b1=b1.split("\t")
            
            xx_1=int(tmp_b1[12].split(":")[-1])
            xy_1=int(tmp_b1[13].split(":")[-1])
            
            for i2,b2 in enumerate(buffer):
                if i1 >= i2:
                    continue
                    
                tmp_b2=b2.split("\t")    
                xx_2=int(tmp_b2[12].split(":")[-1])
                xy_2=int(tmp_b2[13].split(":")[-1])
                
                i_type="I"
                if abs(xy_1-xx_2) <= distance_definition:
                    i_type="D"
                    
                #print(i_type,previous_read_id,tmp_b1[2],tmp_b1[3],tmp_b1[12],tmp_b1[13],tmp_b1[14],tmp_b1[15],tmp_b1[16],tmp_b1[17],tmp_b2[2],tmp_b2[3],tmp_b2[12],tmp_b2[13],tmp_b2[14],tmp_b2[15],tmp_b2[16],tmp_b2[17],sep="\t")
        
        writeMatrix(header_rows,header_cols,matrix,"chimeraTie.matrix.gz")
        
    bam2itx_proc.wait()
    
    verboseprint("")

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
    
    tmp = file.split("\.")
    del tmp[-1]

    file_name = '.'.join(tmp)

    return file_name

def count_lines(file):
    fh = input_wrapper(file)
    
    count = 0
    for _ in fh:
        count += 1
    return count

def iterate_chimera_tie(iter,fastq_file,sam_file,header_file,bowtie2_idx_prefix,genome_name,bowtie_path,min_seq_len,bowtie2_options):

    sam_fh = output_wrapper(sam_file,True,True)

    # open same file
    bowtie_sam_file=genome_name+'___'+str(iter)+'.bowtie.sam'
    bowtie_unaligned_file=genome_name+'___'+str(iter)+'.unaligned.fastq.gz'
    bowtie_sam_fh = output_wrapper(bowtie_sam_file)

    bowtie_cmd = bowtie_path+ ' '+bowtie2_options+' -x '+bowtie2_idx_prefix+' -U '+fastq_file+' -S '+bowtie_sam_file
    bowtie_args = shlex.split(bowtie_cmd)
    bowtie_proc = subprocess.Popen(bowtie_args,
                        stdout=bowtie_sam_fh,
                        stderr=subprocess.PIPE)
    bowtie_proc.wait()
    _,bowtie_results=bowtie_proc.communicate();
    print(bowtie_results)
    bowtie_sam_fh.close()   

    if not debug and iter > 0:
        os.remove(fastq_file)
        
    # open fastq file
    fastq_file=genome_name+'__'+str(iter)+'.fastq.gz'
    fastq_fh = output_wrapper(fastq_file)
    
    # open header file
    header_fh = output_wrapper(header_file,True,True)
    
    # read sam file
    bowtie_sam_fh=input_wrapper(bowtie_sam_file)
    
    for i,line in enumerate(bowtie_sam_fh):
        line=line.rstrip("\n")
        if line.startswith("#") or line.startswith("@"):
            if iter == 0:
                print(line,file=header_fh)
            continue
        
        skipSam=False

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
            values = pattern.split(cigar)[:-1]
            cigar_tup=zip(values[0::2],values[1::2])
            qual=qual[::-1]
            seq=revcomp(seq)
        
        for i in cigar_tup:
            cigar_dict[i[1]]+=int(i[0])
        
        # check later to ensure this is correct
        
        # length of read
        read_length=cigar_dict['M']+cigar_dict['I']+cigar_dict['S']+cigar_dict['=']+cigar_dict['X']
        # match length of cigar
        match_length=cigar_dict['M']+cigar_dict['=']+cigar_dict['X']+cigar_dict['N']+cigar_dict['I']
        # size in ref
        span_length=cigar_dict['M']+cigar_dict['=']+cigar_dict['X']+cigar_dict['N']+cigar_dict['D']
        
        offset=None
        if(len(qname.split(":::")) == 2):
            offset=qname.split(":::")[-1]
        qname=qname.split(":::")[0]
        offset_start=offset_end=0
        if(offset != None):
            offset_start,offset_end=offset.split("-")
  
        xx=int(offset_start)
        xx_rel=0
        
        if(len(cigar_tup) > 1):
            left_cigar=cigar_tup[0]
            right_cigar=cigar_tup[-1]    
            
            if(left_cigar[1] == 'S'):
                left_seq=seq[0:int(left_cigar[0])]
                left_qual=qual[0:int(left_cigar[0])]
                left_start=int(offset_start)
                left_end=int(offset_start)+int(left_cigar[0])
                xx += int(left_cigar[0])
                xx_rel += int(left_cigar[0])

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
        xy_rel=xx_rel+match_length
        middle_seq=seq[xx_rel:xy_rel]
        middle_qual=qual[xx_rel:xy_rel]
        middle_start=xx
        middle_end=xy

        # capture every line to aggregrate sam
        tmp=line.split("\t")
        tmp[0]=tmp[0].split(":::")[0]

        if ('S' in cigar) and (int(tmp[4]) == 0) and (len(middle_seq) > min_seq_len):
            skipSam=True
            print("@"+qname+":::"+str(middle_start)+"-"+str(middle_end),middle_seq,"+",middle_qual,sep="\n",file=fastq_fh)

        if cigar != "*":
            tmp.insert(12,"ZR:i:"+str(read_length))
            tmp.insert(12,"ZS:i:"+str(span_length))
            tmp.insert(12,"ZM:i:"+str(match_length))
            tmp.insert(12,"XQ:i:"+str(iter))
            tmp.insert(12,"XY:i:"+str(xy))
            tmp.insert(12,"XX:i:"+str(xx))

        line='\t'.join(tmp)

        # only write to sam, if iter == 0, or valid alignment, or if skipSam not set (a split read)
        if (iter == 0 or cigar != "*") and not skipSam:
            print(line,file=sam_fh)

    bowtie_sam_fh.close()
    fastq_fh.close()
    header_fh.close()

    if not debug:
        os.remove(bowtie_sam_file)

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
    if outfile.endswith('.fastq.gz'):
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
