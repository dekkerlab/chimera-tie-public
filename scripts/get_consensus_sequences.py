#!/bin/env python3

import argparse
import os
import gzip
import copy

from chimera_lib.gtf import GtfFile

#########################################################################
#
#  TODO output ALL possible exon-exon junctions in a separate file
#  pay attention to the coordinates (the coordinates must be with respect to)
#  the gene coordinates and the genome!      mev42lana



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

#################################################################

def get_consensus_exons(input_exon_list):
    if len(input_exon_list) <= 1:
        return input_exon_list

    unsorted_exonlist = copy.deepcopy(input_exon_list)

    exon_list = sorted(unsorted_exonlist, key = (lambda x : x[0]) )

    exists_overlapping_exons = False
    current_exon_index = -1

    while current_exon_index < len(exon_list) - 1:
        if not exists_overlapping_exons:
            current_exon_index += 1
        current_exon = exon_list[current_exon_index]
        exists_overlapping_exons = False

        for i in range( current_exon_index + 1 , len(exon_list) ):
            other_exon = exon_list[i]
            if other_exon[0] >= current_exon[0] and \
               other_exon[0] <= current_exon[1]:
                exists_overlapping_exons = True
                new_exon = (  current_exon[0] ,
                              max( current_exon[1], other_exon[1] ) )
                exon_list.remove(current_exon)
                exon_list.remove(other_exon)
                exon_list = [new_exon] + exon_list
                exon_list = sorted(exon_list, key = (lambda x : x[0]) )
                break

    return exon_list

#################################################################

def subtract_intervals(minuend, subtrahend):
    '''
    Bot minuend and subtrahend are intervals of the form
    minuend = (x,y)
    subtrahend = (a,b)
    It returns the largest portion of minuend that does not
    intersect with subtrahend
    returns minuend - subtrahend
    If the subtrahend covers the whole minuend then
    empty list is returned

    '''
    x, y = minuend[0],    minuend[1]
    a, b = subtrahend[0], subtrahend[1]

    if( x > y ):
        raise(Exception("Incorrect minuend x > y : " + str(minuend)))
    if( a > b ):
        raise(Exception("Incorrect subtrahend a > b : " + str(subtrahend)))

    # if minuend and subtrahend are disjoint,
    # return the minuend
    if (y < a and y < b) or (x > a and x > b ) :
        return [minuend]

    # if minuend is totally inside the subtrahend
    # return false, meaning that the result is the empty interval
    if a <= x and b >= y:
        return []

    # If subtrahend is inside minuend, return two intervals
    # on both sides
    if x < a and y > b:
        return [ (x, a - 1), (b+1, y) ]

    if x <= a and a <= y and y <= b :
        return [(x, a-1)]

    if a <= x and x <= b and b <= y:
        return [(b + 1, y)]

    print("Warning: This functions shouldn't have printed this."
    " There is a problem with the logical flow!")
    print("minuend:", minuend)
    print("subtrahend", subtrahend)


#################################################################

### TODO
### Test this function!!!

def get_consensus_CDS(consensus_exons, consesnus_UTR):
    previous_exons = consensus_exons
    consensus_cds = consensus_exons

    for utr in consesnus_UTR:
        previous_exons = consensus_cds
        consensus_cds = list()

        for exon in previous_exons:
            consensus_cds += subtract_intervals( exon, utr )

    return consensus_cds

#################################################################

def arrange_gtf_contents(gtf_contents):
    '''
    Input is a dictionary of where each key is a gene name and
    the value is its contents
    This function determines the 5' and 3' consensus UTRS
    Consensus exons and repors CDS exons and UTRS separately
    '''

    for gene, contents in gtf_contents.items():
        consensus_UTRS = get_consensus_exons(contents.UTRS)
        consensus_exons = get_consensus_exons(contents.exons)
        consensus_CDS   = get_consensus_CDS(consensus_exons, consensus_UTR)

        # If we found less than, or more than, 2 UTRS print a warning
        # and do not report any UTRs for that gene
        if len(consensus_UTRS) == 2:
            contents.consensus_UTRS = consensus_UTRS
        else:
            contents.consensus_UTRS = tuple()
            print("Warning: " + str(len(consensus_UTRS)) + \
                  " UTRS have been detected for the gene" + gene + \
                  ". So no actual UTRs have been reported.")
        contents.consensus_exons = consensus_exons
        contents.consensus_CDS_exons = consensus_CDS

    return gtf_contents

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

    gtf_contents = arrange_gtf_contents(gtf_contents)

    return gtf_contents

#############################################################

def main():
    arguments = get_arguments()

    gtf_contents = get_gtf_contents(arguments.ig)

    for key, val in gtf_contents:
        print(key, ":\n")
        print(val, "\n", "--------")

###############################################################

def test_get_consensus_exons(  ):
    print("Testing consensus exons function")
    input_0 = list()
    output_0 = get_consensus_exons(input_0)
    print(input_0, "\n", output_0, "\n----------\n")

    input_1 = [  (3,7), (5,11) , (9,15), (20,40), (50,60) ]
    output_1 = get_consensus_exons(input_1)
    print(input_1, "\n", output_1, "\n----------\n")

    input_2 = [(9,15), (20,40),  (3,7), (5,11) , (50,60) ]
    output_2 = get_consensus_exons(input_2)
    print(input_2, "\n", output_2, "\n----------\n")

    input_3 = [(5,10), (10, 20), (70, 80), (75,85), (90, 100),
               (5,10), (99,200) ]
    output_3 = get_consensus_exons(input_3)
    print(input_3, "\n", output_3, "\n----------\n")

################################################################

def test_subtract_intervals():
    test_pairs = [ ( (1, 3) , (7, 9) ) ,
                   ( (2, 10) , (5, 16) ),
                   ( (30, 100), (50, 80) ),
                   ( ( 30, 40 ) , (20, 90) ),
                   ( (80, 120), (70, 95) ) ]

    for pair in test_pairs:
        print("minuend:", pair[0], "\nsubtrahend:", pair[1])
        print( subtract_intervals(pair[0], pair[1]) )



#################################################################

if __name__ == "__main__":
    #main()
    #test_get_consensus_exons()
    test_subtract_intervals()
