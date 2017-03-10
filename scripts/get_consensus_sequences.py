#!/bin/env python3

import argparse
import os
import gzip

from chimera_lib.gtf import GtfFile

#########################################################################

def get_arguments_helper():
    parser = argparse.ArgumentParser(description=
    '''
    TODO
    Rewrite this part  for this particular script!
    ''')
    parser.add_argument("--ig" ,
                        help = "Input GTF File" ,
                        required = False ,
                        metavar = "Input_GTF_File" ,
                        type = str)
    parser.add_argument("--if" ,
                        help = "Input genome sequence file" ,
                        required = False ,
                        metavar = "input_file_list" ,
                        type = str)
    parser.add_argument("--og" ,
                        help = "Output GTF File" ,
                        required = False ,
                        metavar = "Output_GTF_File" ,
                        type = str)
    parser.add_argument("--of" ,
                        help = "Output consensus sequence File" ,
                        required = False ,
                        metavar = "Output_consensus_sequence_File" ,
                        type = str)


    parser.add_argument("--verbose" ,
                        help = "verbose mode" ,
                        required = False ,
                        dest = 'verbose',
                        action = 'store_true')

    return parser.parse_args()

############################################################

def verbose_print(*args, **kwargs):
    pass

############################################################

def get_arguments():
    arguments = get_arguments_helper()

    return arguments

############################################################

def myopen(file, mode='r'):
    if file.lower()[-3:] == ".gz":
        return gzip.open(file, mode + 't' )
    return open(file, mode)

#############################################################

def get_gtf_contents(input_gtf_file):
    # keys are gene names and values are dictionaries of the following form:
    # chr: chromosome
    # strand: + or -
    # start_position: start_position_of_the_gene
    # end_position: en position of the gene
    # exons: list of exons where each entry is a pair of the form
    #          (exon_start, exon_end)
    #        Note that both coordinates are inclusive and 1-based
    gtf_contents = dict()

    gtf_file = GtfFile(input_gtf_file)

    for entry in gtf_file:
        this_gene = entry.attribute_contents["gene_name"]
        if gtf_contents.get(this_gene) == None:
            # initialize the gene entry
            gtf_contents[this_gene] = { "chr"    : entry.seqname,
                                            "strand" : entry.strand,
                                            "start"  : entry.start,
                                            "end"    : entry.end,
                                            "exons"  : list(),
                                            "UTRS" : list(), }

            # Read all UTRS in the UTRs component
            # At the end, walk over the UTRS and if they overlap merge the overlapping pieces and
            # Remove the pieces from the list and put the merged piece in the list
            # eventually we hould end up with two entries in the UTRS list which
            # will give us 5' UTR and  3'UTRS
            # based on the strand and the position of the UTRS
            # we can decide which one is which
                                      
         if entry.start < gtf_contents[this_gene].start:
             gtf_contents[this_gene].start = entry.start
         if entry.end > gtf_contents[this_gene].end:
            gtf_contents[this_gene].end = entry.end
         if entry.feature == "exon":
             this_exon = (entry.start, entry.end)
             if this_exon not in gtf_contents[this_gene]["exons"]:
                 gtf_contents[this_gene]["exons"].append(this_exon)

    return gtf_contents

#############################################################

def main():
    arguments = get_arguments()

    gtf_contents = get_gtf_contents(arguments.ig)

    for key, val in gtf_contents:
        print(key, ":\n")
        print(val, "\n", "--------")


###############################################################

if __name__ == "__main__":
    main()
