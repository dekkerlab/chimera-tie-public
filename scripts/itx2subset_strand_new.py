#!/usr/bin/env python

"""
Chimera-tie
Hakan Ozadam & Mihir Metkar
"""

from __future__ import print_function
from __future__ import division
import argparse
from collections import defaultdict
from datetime import datetime

from chimera_lib.gff import GFFReader

__version__ = "0.0"

def get_commandline_arguments():
    description = '''
    Given an interaction file, whose columns are tab separated, and a GFF file, this script adds
    gene information at the end of the interaction file for each interacting fragment.

    A line of the interation file is of the form:
    inter_type  read_id  chr_1 pos_1 samflag_1 strand_1  XX:i:22 XY:i:176        XQ:i:1  ZM:i:154        ZS:i:154        ZR:i:176
                         chr_2 pos_2 samflag_2 strand_2  XX:i:22 XY:i:176        XQ:i:1  ZM:i:154        ZS:i:154        ZR:i:176

    '''
    parser = argparse.ArgumentParser(description='Extract direcect/indirect interactions from chimeraTie BAM file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--sorted_itx', dest='sorted_itx_file', type=str, required=True,
                        help='input sorted itx file - output of bam2itx.py')
    parser.add_argument('-g', '--gene_annotation', dest='gene_annotation', type=str,
                        required=True, help='path to gene annoation GFF file')
    parser.add_argument('-s', '--strandness', dest='strand', type=str,
                        required=True, help='strandness of the RNA-Seq Library'
                                            'The fragments are going to be matched to the gff entries'
                                            'according to this parameter. If you do not want to'
                                            'consider strand, set it to N'
                                            'F: Forward, R: Reverse, N: None',
                        choices = ('F', 'R', 'N'))
    parser.add_argument('--debug', dest='debug', action='store_true', help='debug mode')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser.parse_args()

#####################################################################################################

def main():
    command_line_arguments = get_commandline_arguments()
    my_gff_read_file = GFFReader(command_line_arguments.gene_annotation)
    # contents_by_chrom attribute contains named tuples sorted by the start position
    print(my_gff_read_file)

    search_result = my_gff_read_file.search_region_of_fragment(frag_start = 73040490, frag_end = 73040495, frag_chrom = 'X',
                                                               frag_strand = '-', strandness = 'N')

    search_result = my_gff_read_file.search_region_of_fragment(frag_start=3, frag_end=5, frag_chrom='X',
                                                               frag_strand='-', strandness='N')

    print(search_result)


if __name__=="__main__":
    main()
