#!/usr/bin/env python

# encoding: utf-8

"""
chimera-tie.py
Created by Bryan R. Lajoie on 09/16/2015
Edited by Mihir Metkar
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

    parser=argparse.ArgumentParser(description='Extract direcect/indirect interactions from chimeraTie BAM file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--sorted_itx', dest='sorted_itx_file', type=str, required=True, help='input sorted itx file - output of bam2itx.py')
    parser.add_argument('-g', '--gene_annotation', dest='gene_annotation', type=str, required=True, help='path to gene annoation GFF file')
    parser.add_argument('--dd', '--distannce_definition', dest='distance_definition', type=int, default=0)
    parser.add_argument('--debug', dest='debug', action='store_true', help='debug mode')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    sorted_itx_file=args.sorted_itx_file
    gene_annotation=args.gene_annotation
    distance_definition=args.distance_definition
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

    verboseprint("\n",end="")

    # check that bowtie2 and samtools are installed
    check_samtools()

    # load GFF file
    verboseprint("processing GFF file ... ")
    genes,gene_header_file=load_gff(gene_annotation)

    itx_name=get_file_name(sorted_itx_file)
    gff_name=get_file_name(gene_annotation)
    itx_name=itx_name+'__'+gff_name

    verboseprint("ensuring itx file is sorted ... ")

    itx_in_fh=input_wrapper(sorted_itx_file)

    last_chr=last_start=last_strand=None
    num_itx=0
    for line_num,line in enumerate(itx_in_fh):
        x=line.rstrip("\n").split("\t")

        if line.startswith("#"):
            continue
        if line.startswith("@"):
            continue

        chr=x[2]
        start=int(x[3])

        if last_chr != None:
            if chr < last_chr:
                sys.exit('must supply sorted itx file! [run bam2itx.py again!]\n\n['+last_chr+'\t'+last_start+'] ['+chr+'\t'+start+']')
            if chr == last_chr and start < last_start:
                sys.exit('must supply sorted itx file! [run bam2itx.py again!]\n\n['+last_chr+'\t'+last_start+'] ['+chr+'\t'+start+']')

        last_chr=chr
        last_start=start

        num_itx += 1

    itx_in_fh.close()

    verboseprint("\t",num_itx," interactions",sep="")

    verboseprint("")

    pc_resolution=int(math.floor(num_itx/10000))

    col1_overlapped_itx_file=overlapitx(itx_name,sorted_itx_file,genes,1,2,5,3,9)
    n_col1_sorted=count_lines(sorted_itx_file)

    verboseprint("")

    col2_sorted_itx_file=sortitx(itx_name,col1_overlapped_itx_file,2,'-k13,13 -k14,14n') #Is sorting proper? (chr_loc2, start_loc2 of _strand)
    n_col1_overlapped=count_lines(col1_overlapped_itx_file)
    os.remove(col1_overlapped_itx_file)

    col2_overlapped_itx_file=overlapitx(itx_name,col2_sorted_itx_file,genes,2,12,15,13,19)
    n_col2_sorted=count_lines(col2_sorted_itx_file)
    os.remove(col2_sorted_itx_file)

    n_col2_overlapped=count_lines(col2_overlapped_itx_file)

    verboseprint("")

    verboseprint("num col1 sorted =",n_col1_sorted)
    verboseprint("num col1 overlapped =",n_col1_overlapped)
    verboseprint("num col2 sorted =",n_col2_sorted)
    verboseprint("num col2 overlapped =",n_col2_overlapped)

    filenames = [gene_header_file, col2_overlapped_itx_file]
    out_fh=output_wrapper(itx_name+'.itx')
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                out_fh.write(line)
    out_fh.close()

    os.remove(col2_overlapped_itx_file)

    verboseprint("")

def overlapitx(prefix,itx_file,genes,col_num,chr_index=2,strand_index=5,start_index=3,matchlength_index=7):

    verboseprint("overlapping itx [col"+str(col_num)+"] with GFF ... ")

    overlapped_itx_file=prefix+'.col'+str(col_num)+'.overlapped.itx'
    itx_out_fh=output_wrapper(overlapped_itx_file,suppress_comments=True)

    c=0
    for i1,i2 in sweep_overlap(itx_file,genes,chr_index,strand_index,start_index,matchlength_index):
        tmp_gene_list=[i1[0],i1[1]['chrom'],i1[1]['strand'],i1[1]['start'],i1[1]['end']]
        i0=i2+tmp_gene_list
        overlapped_sam_line="\t".join(str(x) for x in i0)
        print(overlapped_sam_line,file=itx_out_fh)
        c+=1

    itx_out_fh.close()

    verboseprint("\tdone")

    return(overlapped_itx_file)

def sortitx(prefix,itx_file,col_num=1,sort_opts="",output_file=None):

    verboseprint("sorting itx [col"+str(col_num)+"] ... ")
    if output_file == None:
        output_file=prefix+'.col'+str(col_num)+'.sorted.itx'

    itx_in_fh=input_wrapper(itx_file)
    itx_out_fh=output_wrapper(output_file,suppress_comments=True)

    sort_cmd = "sort "+sort_opts

    sort_args = shlex.split(sort_cmd)
    sort_proc = subprocess.Popen(sort_args,
                    stdin=itx_in_fh,
                    stdout=itx_out_fh,)
    sort_proc.wait()

    itx_in_fh.close()
    itx_out_fh.close()

    verboseprint("\tdone")

    return output_file

def sweep_overlap(file,genes,chr_index=2,strand_index=5,start_index=3,matchlength_index=9):
    """
    invoke a 'sweeping' algorithm that requires position-sorted file for determing overlap
    """

    itx_fh=open(file,"r")

    get_gene_pos = ( lambda x: (x[1]["chrom"],x[1]["strand"],int(x[1]["start"]),int(x[1]["end"])) )
    get_sam_pos = ( lambda x: (x[chr_index],x[strand_index],int(x[start_index]),int(int(x[start_index])+int(x[matchlength_index].split(":")[-1]))) )

    gene_iter=(i for i in genes)
    itx_iter=(i.rstrip("\n").split("\t") for i in itx_fh)

    c=0
    for i1,i2 in intersection_iter(gene_iter,itx_iter,get_gene_pos,get_sam_pos):
        yield i1,i2
        c=c+1

def intersection_iter(loc1_iter,loc2_iter,posf1,posf2):

    loc2_buffer=[]

    for loc1 in loc1_iter:

        loc1_chr,loc1_strand,loc1_start,loc1_end=posf1(loc1)

        if loc1_start>loc1_end:
            sys.exit('loc1 start>end: '+str((loc1_chr,loc1_start,loc1_end,loc1))+')')

        # remove from buffer locations that have been passed

        new_loc2_buffer=[]

        for i in loc2_buffer:
            if i!=None:
                i_chr,i_strand,i_start,i_end=posf2(i)
            if i==None or i_chr>loc1_chr or (i_chr==loc1_chr and i_end>=loc1_start): #do i need strand info?
                new_loc2_buffer.append(i)

        loc2_buffer=new_loc2_buffer

        # add to buffer locations that intersect

        while True:

            if len(loc2_buffer)>0:

                if loc2_buffer[-1]==None:
                    break

                last_chr,last_strand,last_start,last_end = posf2(loc2_buffer[-1])

                if last_chr>loc1_chr: #do i need strand info?
                    break

                if last_chr==loc1_chr and last_start>loc1_end:
                    break

            try:

                newloc2=loc2_iter.next()

                newloc2_chr,newloc2_strand,newloc2_start,newloc2_end=posf2(newloc2)

                if newloc2_start>newloc2_end:
                    sys.exit('loc2 start>end: '+str((newloc2_chr,newloc2_start,newloc2_end)))

                # add location to buffer if relevant
                if newloc2_chr==None or newloc2_chr>loc1_chr or (newloc2_chr==loc1_chr and newloc2_end>=loc1_start):
                    loc2_buffer.append(newloc2)

            except StopIteration: # if loc2_iter ended

                loc2_buffer.append(None)

        # yield loc1 x loc2_buffer

        for loc2 in loc2_buffer[:-1]:
            yield loc1,loc2



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

def is_overlap(a, b):
    """test to for overlap between two intervals.
    """

    if(a[0] > a[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(a[0])+' > end '+str(a[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    if(b[0] > b[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(b[0])+' > end '+str(b[1])+'!\n\t'+str(a)+' '+str(b)+'\n')

    if a[0] < b[0] and a[1] > b[1]:
        return((b[1]-b[0])+1)

    if b[0] < a[0] and b[1] > a[1]:
        return((a[1]-a[0])+1)

    if b[0] < a[0]:
        a,b=flip_intervals(a,b)

    return max(0, ( min(a[1],b[1]) - max(a[0],b[0]) ) )

def load_gff(gene_annotation):

    if not os.path.isfile(gene_annotation):
        die("GFF file is missing"+gene_annotation)

    gene_file_name=get_file_name(gene_annotation)

    genes=dict()
    ignored_genes=set()
    all_genes=set()

    header2index = dict()

    gff_fh=input_wrapper(gene_annotation)

    n_gff=0
    for i,line in enumerate(gff_fh):
        line=line.rstrip("\n")
        x=line.split("\t")

        if(len(x) == 0):
            continue

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
        strand=x[header2index["strand"]]
        start=int(x[header2index["txStart"]])
        end=int(x[header2index["txEnd"]])
        # build in exon info here - consensus exons
        # to do
        #

        n_gff+=1
        all_genes.add(name2)

        if name2 in ignored_genes:
            continue
        elif name2 not in genes:
            genes[name2]=defaultdict(list)
            genes[name2]["chrom"]=chrom
            genes[name2]["strand"]=strand
            genes[name2]["start"]=start
            genes[name2]["end"]=end
        else:
            if chrom != genes[name2]["chrom"]:
                verboseprint("WARNING: ignoring "+name2+" (exists on multiple chromosomes!)")
                ignored_genes.add(name2)
                del genes[name2]
                continue

            if start < genes[name2]["start"]:
               genes[name2]["start"]=start
            if end > genes[name2]["end"]:
               genes[name2]["end"]=end


    n_genes=len(genes)
    n_ignored_genes=len(ignored_genes)
    n_all_genes=len(all_genes)
    verboseprint("\t",n_gff," entires in supplied GFF",sep="")
    verboseprint("\t",n_all_genes," genes in supplied GFF",sep="")
    verboseprint("\tignored",n_ignored_genes,"multi-chr genes in GFF")
    verboseprint("\tkept",n_genes,"genes in GFF")

    sorted_genes=sorted(genes.items(), key=lambda x: (x[1]['chrom'],int(x[1]['start'])))

    n_overlapped_genes=0
    last_gene=last_gene_name=last_gene_chrom=last_gene_strand=last_gene_start=last_gene_end=None
    for i in sorted_genes:

        if last_gene != None:
            last_gene_name=last_gene[0]
            last_gene_chrom=last_gene[1]["chrom"]
            last_gene_strand=last_gene[1]["strand"]
            last_gene_start=last_gene[1]["start"]
            last_gene_end=last_gene[1]["end"]

        name=i[0]
        chrom=i[1]["chrom"]
        strand=i[1]["strand"]
        start=i[1]["start"]
        end=i[1]["end"]

        overlap=None
        if last_gene_chrom == chrom and last_gene_strand == strand:
            print(name,last_gene_strand,strand,start,end)
            overlap=is_overlap((last_gene_start,last_gene_end),(start,end))

        if overlap > 0:
            n_overlapped_genes+=2
            if last_gene_name in genes:
                del genes[last_gene_name]
            if name in genes:
                del genes[name]

        last_gene=i

    verboseprint("\tremoved",n_overlapped_genes,"overlapping genes!")

    n_genes=len(genes)
    verboseprint("\tkept",n_genes,"genes in GFF")

    genes=sorted(genes.items(), key=lambda x: (x[1]['chrom'],int(x[1]['start'])))

    gene_header_file=gene_file_name+".headers"
    gene_fh=output_wrapper(gene_header_file)

    for i in genes:
        name=i[0]
        chrom=i[1]["chrom"]
        strand=i[1]["strand"]
        start=i[1]["start"]
        end=i[1]["end"]
        length=end-start

        print("@",name,"\t",chrom,"\t",strand,"\t",start,"\t",length,sep="",file=gene_fh)

    gene_fh.close()

    verboseprint("")

    return(genes,gene_header_file)

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

def check_samtools():
    samtools_bin = "samtools"

    global samtools_path
    samtools_version = None
    samtools_path=which(samtools_bin)
    if samtools_path == None:
           die("Error: samtools not found on this system.")

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
