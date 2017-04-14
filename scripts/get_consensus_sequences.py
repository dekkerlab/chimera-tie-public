#!/bin/env python3

import argparse
import os
import gzip
import copy
from time import gmtime, strftime

from chimera_lib.gtf import GtfFile
from chimera_lib.fasta import FastaFile, FastaEntry,\
                              reverse_complement

#########################################################################
#
#  TODO output ALL possible exon-exon junctions in a separate file
#  pay attention to the coordinates (the coordinates must be with respect to)
#  the gene coordinates and the genome!

# TODO
# Organize tests in this script.
# If possible, arrange them in unit tests.

def get_arguments_helper():
    parser = argparse.ArgumentParser(description=
    '''
    TODO
    Rewrite this part  for this particular script!
    ''')
    parser.add_argument("--ig" ,
                        help = "Input GTF File" ,
                        required = True ,
                        metavar = "Input_GTF_File" ,
                        type = str)
    parser.add_argument("--if" ,
                        help = "Input genome sequence file" ,
                        required = True ,
                        metavar = "input_genome_sequence" ,
                        dest = "fasta",
                        type = str)
    parser.add_argument("--of" ,
                        help = "Output consensus sequence File Prefix" ,
                        required = True ,
                        metavar = "Output_consensus_sequence_File_Prefix" ,
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

def pre_verbose_print(*args, **kwargs):
    time_string = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    return print(time_string + " >", *args, **kwargs)

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

    exon_list = sorted( unsorted_exonlist, key = (lambda x : x[0]) )

    exists_overlapping_exons = False
    current_exon_index = -1

    while current_exon_index < (len(exon_list) - 1):
        if not exists_overlapping_exons:
            current_exon_index += 1
        current_exon = exon_list[current_exon_index]
        exists_overlapping_exons = False

        for i in range( current_exon_index + 1 , len(exon_list) ):
            other_exon = exon_list[i]
            # If there are overlapping exons, merge them
            if other_exon[0] >= current_exon[0] and \
               other_exon[0] <= current_exon[1]:
                exists_overlapping_exons = True
                new_relative_start = other_exon[0] - current_exon[0]
                new_relative_end   = other_exon[1] - current_exon[0]
                relative_starts = copy.deepcopy(current_exon[2])
                relative_ends = copy.deepcopy(current_exon[3])
                if new_relative_start not in relative_starts:
                    relative_starts.append(new_relative_start)
                if new_relative_end not in relative_ends:
                    relative_ends.append(new_relative_end)
                new_exon = (  current_exon[0] ,
                              max( current_exon[1], other_exon[1] ),
                              relative_starts, relative_ends )
                exon_list.remove(current_exon)
                exon_list.remove(other_exon)
                exon_list = [new_exon] + exon_list
                # Not very efficient but not too slow
                # We should fix this soon
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
        return [ (x, a - 1, [0] , [ (a - 1) - x]), (b+1, y, [0], [(y-b)-1]) ]

    if x <= a and a <= y and y <= b :
        return [(x, a-1, [0], [(a-1) - x])]

    if a <= x and x <= b and b <= y:
        return [(b + 1, y, [0], [(b + 1) - y])]

    print("Warning: This functions shouldn't have printed this."
    " There is a problem with the logical flow!")
    print("minuend:", minuend)
    print("subtrahend", subtrahend)

#################################################################
def get_consensus_CDS(consensus_exons, consesnus_UTR):
    '''
    ### TODO
    ### Test this function!!!
    ###
    ### TODO
    ### Possible performance improvements
    ### If the exons are ordered we don't need to perform every subtraction
    ### For starters, we are being on the cautious side and
    ### doing the subtraction for each case
    '''
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
    verbose_print("Arranging gtf contents...")

    for gene, contents in gtf_contents.items():
        consensus_UTRS = get_consensus_exons(contents["UTRS"])
        consensus_exons = get_consensus_exons(contents["exons"])

        # If we found less than, or more than, 2 UTRS print a warning
        # and do not report any UTRs for that gene
        if len(consensus_UTRS) == 2:
            contents["consensus_UTRS"] = consensus_UTRS
        else:
            contents["consensus_UTRS"] = tuple()
            print("Warning: " + str(len(consensus_UTRS)) + \
                  " UTRS have been detected for the gene" + gene + \
                  ". So no actual UTRs have been reported.")

        consensus_CDS   = get_consensus_CDS(consensus_exons,
                                            contents["consensus_UTRS"])
        contents["consensus_exons"] = consensus_exons
        contents["consensus_CDS_exons"] = consensus_CDS

    return gtf_contents

#############################################################
def get_sequences_from_consensus_list(gtf_contents, genome_fasta_file,
                                      output_file,
                                      UTR_only = False):
   '''
   Given the coordinates of the consensus exons or CDS_exons and
   the whole genome sequence, this function puts together the sequence of the
   consensus transcripts
   '''


#############################################################
def get_gtf_contents(input_gtf_file):
    # keys are gene names and values are dictionaries of the following form:
    # chr: chromosome
    # strand: + or -
    # start_position: start_position_of_the_gene
    # end_position: en position of the gene
    # exons: list of exons where each entry is a pair of the form
    #          (exon_start, exon_end, relative_starts, relative_ends)
    #        Note that both start and end coordinates are inclusive and 1-based
    # relative_starts and relative ends are
    # both with respect to the start of the exon.
    verbose_print("Getting gtf_contents...")
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
                                            "UTRS" : list(),
                                             }

            # Read all UTRS in the UTRs component
            # At the end, walk over the UTRS and if they overlap merge the overlapping pieces and
            # Remove the pieces from the list and put the merged piece in the list
            # eventually we hould end up with two entries in the UTRS list which
            # will give us 5' UTR and  3'UTRS
            # based on the strand and the position of the UTRS
            # we can decide which one is which
        gene_contents = gtf_contents[this_gene]
        if entry.start < gene_contents["start"]:
             gene_contents["start"] = entry.start
        if entry.end > gene_contents["end"]:
            gene_contents["end"] = entry.end
        if entry.feature.lower() == "exon":
             this_exon = (entry.start, entry.end,
                          [0], [entry.end - entry.start] )
             if this_exon not in gene_contents["exons"]:
                 gene_contents["exons"].append(this_exon)
        if entry.feature.lower() == "utr":
             this_utr = (entry.start, entry.end,
                        [0], [entry.end - entry.start] )
             if this_utr not in gene_contents["UTRS"]:
                 gene_contents["UTRS"].append(this_utr)

    gtf_contents = arrange_gtf_contents(gtf_contents)

    return gtf_contents

#############################################################
def get_genome_sequence(fasta_file):
    fasta_contents = dict()
    with FastaFile(fasta_file) as fasta_input:
        for entry in fasta_input:
            fasta_contents[entry.header] = entry.sequence
    return fasta_contents

#############################################################
def assemble_gene_consensus_sequence( exons, strand, chromosome,
                                      genome_sequence ):

    consensus_sequence = ""
    sequence_pieces    = list()
    chr_sequence = genome_sequence[chromosome]

    for exon in exons:
        sequence_pieces.append( chr_sequence[ (exon[0] - 1) : exon[1] ] )

    consensus_sequence = "".join(sequence_pieces)

    if strand == "-":
        consensus_sequence = reverse_complement(consensus_sequence)

    return consensus_sequence

############################################################
def test_assemble_gene_consensus_sequence():
    exons_1 = ( (4, 6), (12, 14), (19, 20) )
    strand_1 = "+"
    chromosome_1 = "chr2"
    genome_sequence_1 = { "chr10" : "AATT",
                        "chr2" : "AAAGATCCCCCTTATTTTGCAAAAAA" }
    expected_output_1 = "GATTTAGC"
    observed_output_1 = assemble_gene_consensus_sequence( exons = exons_1,
                                          strand = strand_1,
                                          chromosome = chromosome_1 ,
                                          genome_sequence = genome_sequence_1 )

    print("Expected Sequence:", expected_output_1)
    print("Observed Sequence:", observed_output_1, "\n")

    expected_output_2 = "GCTAAATC"
    observed_output_2 = assemble_gene_consensus_sequence( exons = exons_1,
                                          strand = "-",
                                          chromosome = chromosome_1 ,
                                          genome_sequence = genome_sequence_1 )
    print("Expected Sequence:", expected_output_2)
    print("Observed Sequence:", observed_output_2)
#############################################################
def write_consensus_sequences(output_sequence_file,
                              junction_file,
                              genome_sequence_file,
                              gtf_contents, cds_only = False):

    if cds_only:
        verbose_output_target = "exons"
    else:
        verbose_output_target = "cds"

    verbose_print("Writing consensus sequences for",
                   verbose_output_target, "...")

    genome_sequence = get_genome_sequence(genome_sequence_file)
    if cds_only:
        exon_selector = "consensus_CDS_exons"
    else:
        exon_selector = "consensus_exons"

    with myopen(output_sequence_file, "w") as output_stream,\
         myopen(junction_file, "w") as junction_stream:
        for gene, contents in gtf_contents.items():
            consensus_sequence = assemble_gene_consensus_sequence(
                                        exons      = contents[exon_selector],
                                        strand     = contents["strand"],
                                        chromosome = contents["chr"],
                                        genome_sequence = genome_sequence)
            this_fasta_entry = FastaEntry(header=gene,
                                          sequence=consensus_sequence)
            print(this_fasta_entry, file = output_stream)

            gene_length = 0
            for exon in contents[exon_selector]:
                gene_length += exon[1] - exon[0] + 1
            exon_exon_junctions = get_exon_exon_junctions(gene,
                                              gene_length,
                                              contents["start"],
                                              contents["end"],
                                              contents["strand"],
                                              contents[exon_selector])
            print(exon_exon_junctions, file = junction_stream)

#############################################################

def get_exon_exon_junctions(gene_name, gene_length,
                            gene_start, gene_end, strand, exons):
    # the junctions in contents are relative to the exon start.
    # We need to make them relative to the gene start considering the
    # strand of the gene.
    all_junctions = list()
    offset = 0

    for exon in exons:
        starts = list(map( lambda x: x + offset, exon[2]) )
        ends   = list(map( lambda x: x + offset, exon[3]))
        if strand == "-":
            length_minus_one = gene_end - gene_start
            starts = list( map( lambda x: gene_length - x - 1, starts) )
            ends   = list( map( lambda x: gene_length - x - 1, ends) )
        exon_length = (exon[1] - exon[0]) + 1
        offset += exon_length

        # the last nucleotide of an exon can be
        # the first nucleotide of another exon.
        # In that case When the two are merged,
        # we will get the same junction listed twice.
        # In order to avoid that, we keep trac of existing junction positions
        # and ignore a junction position if it already exists.
        existing_junction_positions = list()

        for junction_position in (starts + ends):
            if junction_position in existing_junction_positions:
                continue
            entry_contents = (gene_name, junction_position,
            junction_position + 1, gene_name + "_" + str(junction_position),
            0, strand  )
            entry_contents = list( map(str, entry_contents) )
            all_junctions.append( "\t".join(entry_contents) )

    return "\n".join(all_junctions)

#############################################################

#############################################################
def main():
    arguments = get_arguments()
    print(arguments)
    global verbose_print

    if arguments.verbose:
        verbose_print = pre_verbose_print

    gtf_contents = get_gtf_contents(arguments.ig)
    print(gtf_contents)
    exon_sequence_file = arguments.of + "_consensus_exon_sequence.fa"
    exon_exon_junction_file = arguments.of + "_all_exon_exon_junctions.bed"

    write_consensus_sequences( output_sequence_file = exon_sequence_file,
                               junction_file = exon_exon_junction_file,
                               genome_sequence_file = arguments.fasta ,
                               gtf_contents = gtf_contents,
                               cds_only = False )

    cds_sequence_file = arguments.of + "_consensus_cds_sequence.fa"
    cds_exon_exon_junction_file = arguments.of + "_cds_all_exon_exon_junctions.bed"
    write_consensus_sequences( output_sequence_file = cds_sequence_file,
                               junction_file = cds_exon_exon_junction_file,
                               genome_sequence_file = arguments.fasta,
                               gtf_contents = gtf_contents,
                               cds_only = True )

    return 0

###############################################################
# For the final version, we need to write unit tests instead of
# the following adhoc implementation
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
#### TEST Functions ############################################
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
def test_consensus_CDS():
    exons_1 = ((10, 50), (100, 200), (400, 500) )
    utrs_1 = ( (10, 25) , (451, 500) )
    output_1 = get_consensus_CDS(exons_1, utrs_1)
    expected_output_1 = [ (26,50), (100, 200), (400, 450) ]
    print("Observed_output:", output_1)
    print("Expected output:", expected_output_1)
    print("###############################################")

    exons_2 = ((10, 50), (100, 200), (400, 500) )
    utrs_2 = ( (10, 149) , (451, 500) )
    output_2 = get_consensus_CDS(exons_2, utrs_2)
    expected_output_2 = [ (150, 200), (400, 450) ]
    print("Observed_output:", output_2)
    print("Expected output:", expected_output_2)
    print("###############################################")

#################################################################

if __name__ == "__main__":
    main()
    #test_get_consensus_exons()
    #test_subtract_intervals()
    #test_consensus_CDS()
    #fasta_file = "/chimera-tie/chimera-tie/sample_data/"
    #             "test_consensus/sample_1.fa"
    #test_assemble_gene_consensus_sequence()
