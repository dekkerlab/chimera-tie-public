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

    parser=argparse.ArgumentParser(description='Extract direcect/indirect interactions from chimeraTie BAM file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bam', dest='bam_file', type=str, required=True, help='input bam file - output of chimeraTie.py')
    parser.add_argument('--dd', '--distannce_definition', dest='distance_definition', type=int, default=0)
    parser.add_argument('--dedupe', dest='de_dupe', action='store_true', help='de dupe itx file')
    parser.add_argument('--debug', dest='debug', action='store_true', help='debug mode')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    bam_file=args.bam_file
    distance_definition=args.distance_definition
    de_dupe=args.de_dupe
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

    bam_name=get_file_name(bam_file)

    itx_file=bam_name+'.itx'
    itx_fh=output_wrapper(itx_file,suppress_comments=True)

    dupe_file=bam_name+'.dupes.gz'
    dupe_fh=output_wrapper(dupe_file)

    verboseprint("bam2stats ... [",itx_file,"]")
    bam2itx_cmd = "samtools view "+bam_file
    bam2itx_args = shlex.split(bam2itx_cmd)

    # get num of lines in bam
    num_bam_lines=0
    bam2itx_init = subprocess.Popen(bam2itx_args, bufsize=4096,
                    stdout=subprocess.PIPE)
    with bam2itx_init.stdout:
        for i,_ in enumerate(iter(bam2itx_init.stdout.readline, b'')):
            pass
        num_bam_lines=i+1
    bam2itx_init.wait()

    verboseprint("\t",num_bam_lines," bam entries",sep="")

    verboseprint("")

    pc_resolution=int(math.ceil(num_bam_lines/10000))

    dupes=set()

    verboseprint("bam2itx ... [",itx_file,"]")
    bam2itx_proc = subprocess.Popen(bam2itx_args,bufsize=4096,
                    stdout=subprocess.PIPE)

    previous_read_id=None
    buffer=[]
    with bam2itx_proc.stdout:
        for line_num,line in enumerate(iter(bam2itx_proc.stdout.readline, b'')):
            x=line.split("\t")

            if line_num % pc_resolution == 0:
                pc=(line_num/(num_bam_lines-1))*100
                #verboseprint("\r",""*50,"\r\t"+str(line_num)+" / "+str(num_bam_lines-1)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                if verbose: sys.stdout.flush()

            #unmapped or secondary
            if( (int(x[1]) & 0x4) or (int(x[1]) & 0x100) ):
                continue

            current_read_id=x[0]

            if((current_read_id != previous_read_id) and (previous_read_id != None) and (len(buffer) != 0)):

                key_list=[]
                for i,b in enumerate(buffer):
                    tmp_b=b.split("\t")
                    tmp_key=tmp_b[2]+"_"+tmp_b[3]+"_"+tmp_b[12]
                    key_list.append(tmp_key)

                key=":".join(key_list)

                if de_dupe:
                    if key in dupes:
                        print(key,file=dupe_fh)
                        previous_read_id=current_read_id
                        buffer=[]
                    else:
                        dupes.add(key)

                for i1,b1 in enumerate(buffer):
                    tmp_b1=b1.split("\t")
                    strand_info_1=int(tmp_b1[1])
                    if strand_info_1==0:
                        #ZZ method, reads are anti-sense to the RNA
                        strand_1='-'
                    elif strand_info_1==16:
                        #ZZ method, reads are anti-sense to the RNA
                        strand_1='+'
                    #print(strand_info_1,strand_1)

                    xx_1=int(tmp_b1[12].split(":")[-1])
                    xy_1=int(tmp_b1[13].split(":")[-1])

                    for i2,b2 in enumerate(buffer):
                        if i1 >= i2:
                            continue

                        tmp_b2=b2.split("\t")
                        strand_info_2=int(tmp_b2[1])
                        if strand_info_2==0:
                            #ZZ method, reads are anti-sense to the RNA
                            strand_2='-'
                        elif strand_info_2==16:
                            #ZZ method, reads are anti-sense to the RNA
                            strand_2='+'
                        #print(strand_info_2,strand_2)

                        xx_2=int(tmp_b2[12].split(":")[-1])
                        xy_2=int(tmp_b2[13].split(":")[-1])

                        i_type="I"
                        if abs(xy_1-xx_2) <= distance_definition:
                            i_type="D"

                        print(i_type,previous_read_id,tmp_b1[2],tmp_b1[3],strand_info_1,strand_1,tmp_b1[12],tmp_b1[13],tmp_b1[14],tmp_b1[15],tmp_b1[16],tmp_b1[17],tmp_b2[2],tmp_b2[3],strand_info_2,strand_2,tmp_b2[12],tmp_b2[13],tmp_b2[14],tmp_b2[15],tmp_b2[16],tmp_b2[17],sep="\t",file=itx_fh)

                # automatically add singletons
                if(len(buffer) == 1):
                    i_type="S"
                    tmp_b=buffer[0].split("\t")
                    strand_info=int(tmp_b[1])
                    if strand_info==0:
                        strand='+'
                    elif strand_info==16:
                        strand='-'
                    print(i_type,previous_read_id,tmp_b[2],tmp_b[3],strand_info,strand,tmp_b[12],tmp_b[13],tmp_b[14],tmp_b[15],tmp_b[16],tmp_b[17],tmp_b[2],tmp_b[3],strand_info,strand,tmp_b[12],tmp_b[13],tmp_b[14],tmp_b[15],tmp_b[16],tmp_b[17],sep="\t",file=itx_fh)

                buffer=[]

            # add current line to buffer
            buffer.append(line)

            # set previous = current
            previous_read_id=current_read_id

        key_list=[]
        for i,b in enumerate(buffer):
            tmp_b=b.split("\t")
            tmp_key=tmp_b[2]+"_"+tmp_b[3]+"_"+tmp_b[12]
            key_list.append(tmp_key)

        key=":".join(key_list)

        if de_dupe:
            if key in dupes:
                print(key,file=dupe_fh)
                previous_read_id=current_read_id
                buffer=[]
                #continue
            else:
                dupes.add(key)

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

                print(i_type,previous_read_id,tmp_b1[2],tmp_b1[3],tmp_b1[12],tmp_b1[13],tmp_b1[14],tmp_b1[15],tmp_b1[16],tmp_b1[17],tmp_b2[2],tmp_b2[3],tmp_b2[12],tmp_b2[13],tmp_b2[14],tmp_b2[15],tmp_b2[16],tmp_b2[17],sep="\t",file=itx_fh)

    bam2itx_proc.wait()
    itx_fh.close()
    dupe_fh.close()

    verboseprint("")
    verboseprint("")

    output_file=bam_name+'.sorted.itx'
    output_file=sortitx(bam_name,itx_file,1,'-k3,3 -k4,4n',output_file)

    verboseprint("")

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

def sweep_overlap(file,genes,chr_index=2,start_index=3,matchlength_index=7):
    """
    invoke a 'sweeping' algorithm that requires position-sorted file for determing overlap
    """

    itx_fh=open(file,"r")

    get_gene_pos = ( lambda x: (x[1]["chrom"],int(x[1]["start"]),int(x[1]["end"])) )
    get_sam_pos = ( lambda x: (x[chr_index],int(x[start_index]),int(int(x[start_index])+int(x[matchlength_index].split(":")[-1]))) )

    gene_iter=(i for i in genes)
    itx_iter=(i.rstrip("\n").split("\t") for i in itx_fh)

    c=0
    for i1,i2 in intersection_iter(gene_iter,itx_iter,get_gene_pos,get_sam_pos):
        yield i1,i2
        c=c+1

def intersection_iter(loc1_iter,loc2_iter,posf1,posf2):

    loc2_buffer=[]

    for loc1 in loc1_iter:

        loc1_chr,loc1_start,loc1_end=posf1(loc1)

        if loc1_start>loc1_end:
            sys.exit('loc1 start>end: '+str((loc1_chr,loc1_start,loc1_end,loc1))+')')

        # remove from buffer locations that have been passed

        new_loc2_buffer=[]

        for i in loc2_buffer:
            if i!=None:
                i_chr,i_start,i_end=posf2(i)
            if i==None or i_chr>loc1_chr or (i_chr==loc1_chr and i_end>=loc1_start):
                new_loc2_buffer.append(i)

        loc2_buffer=new_loc2_buffer

        # add to buffer locations that intersect

        while True:

            if len(loc2_buffer)>0:

                if loc2_buffer[-1]==None:
                    break

                last_chr,last_start,last_end = posf2(loc2_buffer[-1])

                if last_chr>loc1_chr:
                    break

                if last_chr==loc1_chr and last_start>loc1_end:
                    break

            try:

                newloc2=loc2_iter.next()

                newloc2_chr,newloc2_start,newloc2_end=posf2(newloc2)

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
