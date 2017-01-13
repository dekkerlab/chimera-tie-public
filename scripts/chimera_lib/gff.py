#
#
"""
Chimera-tie
Hakan Ozadam & Mihir Metkar
"""

# This scirpt reads, parses and sorts a GFF file whose lines are in the form:
# #bin	chrom	processed	entry_type	txStart	txEnd	strand	name	exonic_part_number	name2
# 691260	X	dexseq_prepare_annotation.py	aggregate_gene	73012040	73049066	+	ENSG00000270641	0	ENSG00000270641_0

from collections import defaultdict, namedtuple

class GFFReader:
    '''This object is a container for the GFF reference. It is initiated by the gff file.
    It parses the gff file and reads it into a dictionary called contents_by_chromosome
     It also has a search feature function. Given a read fragment, search_region_of_fragment
     function can tell us if the given fragment matches a gene'''

    def __init__(self, file):
        file_contents_by_gene = defaultdict()
        column_names = ['bin', 'chrom', 'processed',
                        'entry_type', 'start', 'end',
                        'strand', 'name_1', 'number', 'name_2']
        gff_entry = namedtuple('gff_entry', column_names)

        self.contents_by_chrom = defaultdict(list)
        start_label_index = column_names.index("start")
        end_label_index = column_names.index("end")

        #Read the file contents into a dictionary whose keys are chromosomes
        line_number = 0
        with open(file, 'r') as input_file:
            for line in input_file:
                line_number += 1

                if line[0] == '#':
                    continue
                line_contents = list(line.rstrip().split())
                if len(line_contents) < len(column_names):
                    print("Warning: File: {file}: Line {line_number}:\n  This line contains too few entries:\n {line} \n".format(
                                file = file, line_number = line_number, line = line ) )
                    continue


                line_contents[start_label_index] = int(line_contents[start_label_index])
                line_contents[end_label_index] = int(line_contents[end_label_index])
                this_gff_entry = gff_entry(*line_contents)
                self.contents_by_chrom[this_gff_entry.chrom].append(this_gff_entry)

        self.__sort_by_position__()
        self.chromosomes = self.contents_by_chrom.keys()

        # For strand specific searches, it makes more sense to separate the gff entries according
        # to their strand
        self.contents_by_chrom_plus_strand  = defaultdict(list)
        self.contents_by_chrom_minus_strand = defaultdict(list)
        for chrom in self.chromosomes:
            self.contents_by_chrom_plus_strand[chrom] = \
                list( filter(lambda x: x.strand == "+", self.contents_by_chrom[chrom]) )
            self.contents_by_chrom_minus_strand[chrom] = \
                list(filter(lambda x: x.strand == "-", self.contents_by_chrom[chrom]) )

    def __sort_by_position__(self):
        '''Sort the self.contents_by_name_1 by the start value in increasing order '''
        for name_1 in self.contents_by_chrom.keys():
            self.contents_by_chrom[name_1].sort(key = (lambda x: x.start))

    def __str__(self):
        # Defined for debug purposes mostly. It dumps the dictionary in a tab separated
        # line-by-line form
        result_list = list()
        for key, name_entries in self.contents_by_chrom.items():
            for individual_gff_entry in name_entries:
                result_list.append("\t".join( map(str, individual_gff_entry) )  )

        return("\n".join(result_list))

    def search_region_of_fragment(self, frag_start, frag_end, frag_chrom, frag_strand, strandness = 'N'):
        '''Implements binary search to look up the corresponding
        region of a given fragment in the gff file'''

        if frag_chrom not in self.chromosomes:
            return None



        chrom_contents = self.contents_by_chrom[frag_chrom]
        #DEBUG
        #print("chrom contents: ", chrom_contents)
        if (strandness == 'F' and frag_strand == "+") or\
            (strandness== 'R' and frag_strand == "-"):
            chrom_contents = self.contents_by_chrom_plus_strand[frag_chrom]
        elif strandness != 'N':
            chrom_contents = self.contents_by_chrom_minus_strand[frag_chrom]

        start_index , end_index = 0, len( chrom_contents )

        if len(chrom_contents) < 1:
            return None

        if strandness == 'F':
            if frag_strand !=  chrom_contents[0].strand:
                return None
        if strandness == 'R':
            if frag_strand ==  chrom_contents[0].strand:
                return None

        # DEBUG
        #print("func parameters: ", frag_start, frag_end, frag_chrom, frag_strand, strandness)
        #print("chrom contents: ", chrom_contents)

        while True:
            middle_index = (start_index + end_index) // 2
            current_guess = chrom_contents[middle_index]
            if current_guess.start <= frag_start and\
                current_guess.end >= frag_end:
                    return current_guess
            elif current_guess.start > frag_start:
                end_index = middle_index
            elif current_guess.end < frag_end:
                start_index = middle_index
            if middle_index == start_index or middle_index == end_index:
                break
        return None
