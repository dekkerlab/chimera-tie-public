#!/usr/bin/env python

"""
bamheader2gff.py
Created by Mihir Metkar 11/08/16
"""

from __future__ import print_function
import argparse

__version__ = "1.0"

"""Eg. Command- python Exp_comp.py -i header.sam
@HD    	VN:1.0 	SO:unsorted
@SQ    	SN:ENSG00000094661|ENST00000209540     	LN:1242
@SQ    	SN:ENSG00000088356|ENST00000202017     	LN:1998

output gff file format
#bin   	name   	chrom  	txStart	txEnd  	name2
1      	gene_name 	gene_chrom 	1      	122    	gene_symbol
"""

def main():

    parser=argparse.ArgumentParser(description='Convert bam headers to gff \
                                                format for chimera-tie (without strand). \
                                                Used when mapping to transcriptome.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bam_headers', dest='bam_headers', type=str, required=True, help='output of: samtools view -H in.bam > header.sam')
    parser.add_argument('-o', '--outfile_name', dest='outfile_name', type=str, required=True, help='name the output file')
    parser.add_argument('-s', '--strand', dest='strand', type=str, required=True, help='name the output file', choices = ('T', 'F'))
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    bam_headers=args.bam_headers
    outfile_name=args.outfile_name
    strand=args.strand

    bam_headers_fh=open(bam_headers,'r')
    gff_fh=open(outfile_name,'w')

    if strand=='F':
        print ('#bin','name','chrom','txStart','txEnd','name2', sep="\t",file= gff_fh)

    else:
        #bin	chrom	processed	entry_type	txStart	txEnd	strand	name	exonic_part_number	name2
        print ('#bin','chrom','processed','entry_type','txStart','txEnd','strand','name','exonic_part_number','name2', sep="\t",file= gff_fh)


    line_count=1
    for i,line in enumerate(bam_headers_fh):
        if not line.startswith("@SQ"):
            continue

        line=line.rstrip("\n")
        line_element=line.split("\t")

        name=name2=chrom=line_element[1].split(":")[-1]
        length=line_element[2].split(":")[-1]
        start=1
        end=1+int(length)
        bin_num=line_count
        processed='bamheader2gff'
        entry_type='transcriptome'
        strand_info="+"
        exonic_part_number='None'

        if strand=='F':
            #bin   	name   	chrom  	txStart	txEnd  	name2
            print(bin_num, name, chrom, start, end, name2,sep="\t",file= gff_fh)
            line_count+=1

        else:
            #bin	chrom	processed	entry_type	txStart	txEnd	strand	name	exonic_part_number	name2
            print (bin_num, chrom, processed, entry_type, start, end, strand_info, name, exonic_part_number, name2, sep="\t",file= gff_fh)
            line_count+=1

    bam_headers_fh.close()
    gff_fh.close()




if __name__=="__main__":
        main()
