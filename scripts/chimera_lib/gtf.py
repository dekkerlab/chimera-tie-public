#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
#
#        Dekker Labaratory
#        Department of Biochemistry and Molecular Pharmacology
#        University of Massachusetts Medical School
#        368 Plantation Street, Albert Sherman Center
#        Worcester, MA 01605, USA
#

##########################################################
from sys import stdin
import re
##########################################################

###########################################################
class GtfEntry:
    '''
    Parses a given array into a GTF file entries.
    '''

    _fields = ["seqname" , "source", "feature" , "start" , "end" ,
               "score" , "strand", "frame" , "attribute"]

    def __init__(self, entry_tuple):

        if len(entry_tuple) < len(self._fields):
            raise(AttributeError("The input tuple of GtfEntry must have",
                                 len(self._fields), "entries." ))

        for key , val in zip(self._fields, entry_tuple):
            setattr(self, key, val)

        self.attribute  = self.attribute.strip()
        self.start , self.end = int(self.start) , int(self.end)

        attribute_contents = self.attribute.rstrip(";").split(';')

        self.attribute_contents = dict()

        for a in attribute_contents:

            contents = a.split()
            if len(contents) < 2:
                continue
            key      = contents[0]
            val      = contents[1]
            self.attribute_contents[key] = val.replace("\"", "")


    ################################################
    def __str__(self):
        line_contents = list()
        for f in self._fields:
            line_contents.append(getattr(self, f))
        return "\t".join( map(str ,line_contents) )


############################################################

class GtfFile:
    '''
    Open a file and reads its lines into a GTF entry object.
    This class implements an iterator for the GTF file so that it can be read
    sequentially like an ordinary file
    '''
    def __init__(self , file):
        if(file):
            self.f = open(file , "r")
        else:
            self.f = stdin

    def __getitem__(self, index):

        for raw_line in self.f:
            line = raw_line.strip()
            if ( line == '' or line[0] == "#"):
                continue
            return GtfEntry(line.strip().split("\t"))

        raise IndexError

    def __del__(self):
        if hasattr(self, 'f') and self.f:
            self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


################################################################

def check_gtf_validity(gtf_file_path):
    '''
    Check if the given gtf file is valid.
    :param gtf_file_path: The path of the gtf file to be examined
    :return: Returns emptry string if the given gtf file is valid. If not returns the first error
    '''

    sample_valid_gtf_line = "chr1	unknown	exon	14970	15038	.	-	.	gene_id \"WASH7P\"; " \
                            "gene_name \"WASH7P\"; transcript_id \"NR_024540\"; tss_id \"TSS7514\";"

    error_string = ""
    number_of_exons = 0
    counter = 0
    with GtfFile(gtf_file_path) as gtf_input_stream:
        for entry in gtf_input_stream:
            if not entry.attribute_contents.get("gene_id"):
                error_string += "----\nError with entry:" +  str(entry) + \
                               "\nThis entry has no \"gene_id\" in its attributes.\n------\n"
                return error_string
            counter += 1
            if entry.feature == "exon":
                number_of_exons += 1

    if number_of_exons < 2:
        error_string += "This file has less than 2 exons!"

    return error_string
