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
import subprocess
import os

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
    parser = argparse.ArgumentParser(description='Extract direct/indirect interactions from chimeraTie BAM file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--itx', dest='itx_file', type=str, required=True,
                        help='input ,sorted, itx file - output of bam2itx.py')
    parser.add_argument('-o', '--annotated_itx', dest='output_itx_file', type=str, required=True,
                        help=' itx file with annotation. Gene names ,start , end, strand added')
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

OFFSET = 10
FRAG_1_START_INDEX = 3
FRAG_1_CHR_INDEX   = 2
FRAG_1_MATCH_LENGTH_INDEX = 9
FRAG_1_STRAND_INDEX = 5



def get_fragment_contents(line_contents, offset ,strandness):
    '''Given the columns in array from an individual line of the interaction file
    get the contents in dictionary form'''
    result = {}
    result["frag_start"] =  int(line_contents[FRAG_1_START_INDEX + offset])
    raw_fragment_match = line_contents[FRAG_1_MATCH_LENGTH_INDEX + offset]
    frag_match_contents = raw_fragment_match.split(":")
    if len(frag_match_contents) < 3:
        raise(IOError("\t".join(line_contents) + ": Problem in match contents field"))
    result["frag_end"]   = int(frag_match_contents[2]) + result["frag_start"]
    result["frag_chrom"] = line_contents[FRAG_1_CHR_INDEX + offset]
    result["frag_strand"] = line_contents[FRAG_1_STRAND_INDEX  + offset]
    result["strandness"]  = strandness

    return result

##############################################################################################

def prepare_output(piece_1, piece_2):

    str_array = map(str, (piece_1.name_2, piece_1.chrom, piece_1.strand, piece_1.start, piece_1.end, \
             piece_2.name_2, piece_2.chrom, piece_2.strand, piece_2.start, piece_2.end))
    return "\t".join(str_array)

##############################################################################################

def get_output_header(gff_reader):
    header_list = list()
    sorted_chromosomes = sorted(gff_reader.contents_by_chrom.keys() )

    for key in sorted_chromosomes:
        value = gff_reader.contents_by_chrom[key]
        key_contents = dict()

        for entry in value:
            this_length =  int(entry.end) - int(entry.start)
            this_line_list = ( "@" + entry.name_2, entry.chrom,
                                      entry.strand, entry.start,
                                      this_length )

            key_contents[entry.name_2] = this_line_list

        sorted_contents = sorted(key_contents.values(), key = (lambda x : x[3])  )
        #print("sorted contents:", sorted_contents)
        sorted_contents_str = list()
        for item in sorted_contents:
            item_str = "\t".join( map(str, item) )
            sorted_contents_str.append( item_str )


        header_list += sorted_contents_str


    return "\n".join(header_list)



##############################################################################################

def main():
    command_line_arguments = get_commandline_arguments()
    my_gff_read_file = GFFReader(command_line_arguments.gene_annotation)
    # contents_by_chrom attribute contains named tuples sorted by the start position

    output_header = get_output_header(my_gff_read_file)

    pre_output_file = command_line_arguments.output_itx_file + ".pre"


    with open(command_line_arguments.itx_file) as input_stream,\
         open(pre_output_file, 'w') as output_stream:

        for line in input_stream:
            line_contents = line.strip().split()
            if len(line_contents) < 12:
                continue
            frag1_contents = get_fragment_contents(line_contents = line_contents,
                                                   offset = 0 , strandness = command_line_arguments.strand)
            #search for the first fragment in the gff file
            annotation_search_result_1 = my_gff_read_file.search_region_of_fragment(**frag1_contents)
            if annotation_search_result_1 is None:
                continue

            frag2_contents = get_fragment_contents(line_contents=line_contents,
                                                   offset=OFFSET, strandness=command_line_arguments.strand)
            # search for the second fragment in the gff file
            annotation_search_result_2 = my_gff_read_file.search_region_of_fragment(**frag2_contents)
            if annotation_search_result_2 is None:
                continue

            output_line = line.strip() + "\t" + prepare_output(annotation_search_result_1, annotation_search_result_2)
            print(output_line, file = output_stream)


    print("Sorting the output...")

    with open(command_line_arguments.output_itx_file, "w") as output_stream:
        print(output_header, file = output_stream)

    with open(command_line_arguments.output_itx_file, "a") as output_stream:

        my_process = subprocess.Popen( [ "sort", \
                         "-k2,2", "-k7,7", "-k17,17", "-V",\
                        pre_output_file  ], stdout = output_stream  )
        sort_out, sort_err = my_process.communicate()
        #print(output_header, file = output_stream)
        print(sort_out, file = output_stream)

    os.remove(pre_output_file)
    print("Done!")




if __name__=="__main__":
    main()




####################################################################################################################
# Test code for main, can be deleted later
# print(my_gff_read_file)
#
# search_result = my_gff_read_file.search_region_of_fragment(frag_start=73040490, frag_end=73040495, frag_chrom='X',
#                                                            frag_strand='-', strandness='N')
#
# print("Must find", search_result)
# search_result = my_gff_read_file.search_region_of_fragment(frag_start=3, frag_end=5, frag_chrom='X',
#                                                            frag_strand='-', strandness='N')
#
# print("Must be none: ", search_result)
