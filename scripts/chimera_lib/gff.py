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
    def __init__(self, file):
        file_contents_by_gene = defaultdict()
        column_names = ['bin', 'chrom', 'processed',
                        'entry_type', 'start', 'end',
                        'strand', 'name_1', 'number', 'name_2']
        gff_entry = namedtuple('gff_entry', column_names)

        self.contents_by_chrom = defaultdict(list)
        start_label_index = column_names.index("start")
        end_label_index = column_names.index("end")

        with open(file, 'r') as input_file:
            for line in input_file:
                if line[0] == '#':
                    continue
                line_contents = list(line.rstrip().split())
                if len(line_contents) < len(column_names):
                    print("Warning: This line contains too few reads:\n", line, '\n' )
                    continue


                line_contents[start_label_index] = int(line_contents[start_label_index])
                line_contents[end_label_index] = int(line_contents[end_label_index])
                this_gff_entry = gff_entry(*line_contents)
                self.contents_by_chrom[this_gff_entry.chrom].append(this_gff_entry)

        self.__sort_by_position__()
        self.chromosomes = self.contents_by_chrom.keys()

    def __sort_by_position__(self):
        '''Sort the self.contents_by_name_1 by the start value in increasing order '''
        for name_1 in self.contents_by_chrom.keys():
            self.contents_by_chrom[name_1].sort(key = (lambda x: x.start))

    def __str__(self):
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
        # Make this part in init func and make it faster!!!
        if strandness == 'F':
            chrom_contents = filter( self.contents_by_chrom[frag_chrom] , lambda x:  x.strand == frag_strand )
        elif strandness == 'R':
            chrom_contents = filter(self.contents_by_chrom[frag_chrom], lambda x: x.strand != frag_strand )

        start_index , end_index = 0, len( chrom_contents )

        while True:
            middle_index = (start_index + end_index) // 2
            print(middle_index)
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





