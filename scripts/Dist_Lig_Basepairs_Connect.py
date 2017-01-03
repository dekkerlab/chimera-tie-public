#!/usr/bin/env python

# encoding: utf-8

"""
DistvLigvBasePair.py
Created by Mihir Metkar 12/04/16
"""

from __future__ import print_function
import argparse
from collections import Counter
import re

def main():

    parser=argparse.ArgumentParser(description='Create XY scatter- 3D distance (PDB) vs ligation junction', \
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-j', '--lig_junc', dest='lig_junc', type=str, required=True, \
                                    help='input ligation junction file. Input of itx2scatterplot_strand.py')
    parser.add_argument('-d', '--dist_pdb', dest='dist_pdb', type=str, \
                                    help='input 3D distnce from pdb file. input of Pymol_distance.py')
    parser.add_argument('-c', '--connect', dest='connect', type=str, \
                                    help='Base-pair information if available in connect (.ct) format') #no "required =" field because its optional
    parser.add_argument('-f', '--fasta', dest='fasta', type=str, \
                                    help='RNA sequence in fasta format')
    parser.add_argument('-o', '--out_file', dest='out_file', type=str, required=True, help='name of output file')


    args=parser.parse_args()

    lig_junc_file = args.lig_junc
    dist_pdb_file = args.dist_pdb
    connect_file = args.connect
    fasta_file = args.fasta
    out_file = args.out_file

    lig_junc_fh = open(lig_junc_file,'r')

    junctions_list, junction_nuc = getjunctions(lig_junc_fh)

    #print (junctions_list)
    #print ("\n")
    #print (junction_nuc)

    if dist_pdb_file:
        distvlig(dist_pdb_file, lig_junc_fh, out_file, junctions_list)
        dist2cworld(dist_pdb_file, out_file)

    if connect_file:
        connect2helix(connect_file, out_file)
        junction_frequency(connect_file, junction_nuc, out_file )

    if fasta_file:
        print ("\nProcessing ligation junctions file and connect file to extract corresponding base pairs...\n")
        lig2other(junctions_list, fasta_file, out_file)

    lig_junc_fh.close()


def getjunctions(lig_junc_file):
    junc_nuc = []
    junctions = []
    for i,each_line in enumerate(lig_junc_file):
        junc = []
        line = each_line.rstrip("\n").split("\t")

        species1 = line[0]
        species2 = line[1]
        nucleotide5 = int(line[2])
        nucleotide3 = int(line[3])
        lig_count = int(line[4])

        junc.append(species1)
        junc.append(species2)
        junc.append(nucleotide5)
        junc.append(nucleotide3)
        junc.append(lig_count)

        junctions.append(junc)

        junc_nuc.append(nucleotide5)
        junc_nuc.append(nucleotide3)
    return (junctions, junc_nuc)


def distvlig(dist_file, lig_junc_file, output_file, junc_list):
    print ("Processing ligation junctions file to extract corresponding distances...")
    print ("\tCan only process intra-RNA ligations...\n")

    dist_pdb_fh = open(dist_file,'r')
    out_lig_dist_fh = open(output_file+'_lig_dist.txt','w')

    dist_list = []
    for i,each_dist in enumerate(dist_pdb_fh):
        if each_dist.startswith("resi"):
            continue
        dist = each_dist.rstrip("\n").split("\t")
        dist_list.append(float(dist[4]))

    min_dist = min(dist_list)
    max_dist = max(dist_list)
    dist_pdb_fh.seek(0)

    print ("\nWrote ligation vs Distance file!\n")
    #The second for loop reads the whole ValuesFile. after its execution, the file pointer is at the end of the file,
    #and there is no more values to be read from it. You should reset the ValuesFile file pointer before starting
    #the second for loop to start reading ValuesFile from the beginning again

    print ("#minimum distance:",min_dist, ", maximum distance:",max_dist, file = out_lig_dist_fh)

    print ("Species1","Species2","resi1","resi2","Distance in 3D (Ang)","Ligation Frequency", sep = "\t", file = out_lig_dist_fh)

    for junct in junc_list:
        species1 = junct[0]
        species2 = junct[1]
        nucleotide5 = int(junct[2])
        nucleotide3 = int(junct[3])
        lig_count = int(junct[4])

        #print (nucleotide5,nucleotide3,lig_count, sep = "\t")

        for i,each_dist in enumerate(dist_pdb_fh):
            if each_dist.startswith("resi"):
                continue
            dist = each_dist.rstrip("\n").split("\t")

            resi1 = int(dist[0])
            resi2 = int(dist[1])
            distance = float(dist[4])

            if nucleotide5 == resi1 and nucleotide3 == resi2:
                print (species1, species2, resi1, resi2, distance, lig_count, sep="\t", file = out_lig_dist_fh)

        dist_pdb_fh.seek(0)

    out_lig_dist_fh.close()

def dist2cworld(dist_file, output_file):
    print ("Processing distances file to convert to pairwise distances (input for cworld script)...")
    print ("\tCan only process intra-RNA ligations...\n")

    dist_pdb_fh = open(dist_file,'r')
    out_dist2pw_fh = open(output_file+'_dist_pairwise.txt','w')

    print ("yHeader", "xHeader", "cScore", sep = "\t", file = out_dist2pw_fh)

    #yHeader xHeader cScore
    #RNA5S1__0|hg19|RNA5S1:1-2       RNA5S1__114|hg19|RNA5S1:115-116 1.00000000

    for i,each_dist in enumerate(dist_pdb_fh):
        if each_dist.startswith("resi"):
            header = each_dist.rstrip("\n").split("\t")
            header_X = header[0]
            header_Y = header[1]

        #Since PDB file doesnt give distances for
            print (header_Y+"__0|hg19|"+header_Y+":1-2", header_X+"__0|hg19|"+header_X+":1-2", "0.0", sep = "\t", file = out_dist2pw_fh)

        else:
            dist = each_dist.rstrip("\n").split("\t")
            X = int(dist[0])
            X_head = header_X+'__'+str((X-1))+'|hg19|'+header_X+':'+str(X)+'-'+str((X+1))
            Y = int(dist[1])
            Y_head = header_Y+'__'+str((Y-1))+'|hg19|'+header_Y+':'+str(Y)+'-'+str((Y+1))
            score = float(dist[4])
            print (Y_head, X_head, score, sep = "\t", file = out_dist2pw_fh)

def connect2helix(input_file, output_file):
    connect_fh = open(input_file,'r')
    connect2helix_fh = open(output_file+'_ct2helix.txt','w')
    first_line = connect_fh.readline().rstrip("\n")
    print ("#"+first_line, file = connect2helix_fh)
    print ("i", "j", "length", "value", sep = "\t", file = connect2helix_fh)

    for line in connect_fh:
        line = line.rstrip("\n").split(" ")
        if int(line[4]) == 0:
            continue
        else:
            print (line[0], line[4], 1, 1, sep = "\t", file = connect2helix_fh)
    connect_fh.seek(0)


def lig2other(junctions_list, fasta_file, out_file):
    fasta_fh = open(fasta_file,'r')
    lig_connect_fh = open(out_file+'_lig2ct.txt','w')
    lig_helix_fh = open(out_file+'_lig2helix.txt','w')
    lig_connect_filtered_fh = open(out_file+'_lig2ct_filtered.txt','w')
    lig_helix_filtered_fh = open(out_file+'_lig2helix_filtered.txt','w')



    #Write a .ct file with ligation junctions instead of base pairs.
    print ("\nProcessing fasta file and ligations junctions file to make a .ct and filtered.ct file...\n")
    print ("\nProcessing fasta file and ligations junctions file to make a .helix and filtered.helix file...\n")

    letters = []
    nuc_num = 0
    for lines in fasta_fh:
        if lines.startswith(">"):
            continue
        lines = lines.rstrip("\n")
        for i,letter in enumerate(lines):
            nuc_num = nuc_num + 1
    max_num_num = nuc_num
    print (max_num_num, sep ="\t", file = lig_connect_fh)
    print (max_num_num, sep ="\t", file = lig_connect_filtered_fh)
    print ('#'+ str(max_num_num), file = lig_helix_fh)
    print ('#'+ str(max_num_num), file = lig_helix_filtered_fh)
    print ("i","j", "length", "value", sep = "\t", file = lig_helix_fh)
    print ("i","j", "length", "value", sep = "\t", file = lig_helix_filtered_fh)
    fasta_fh.seek(0)

    letters = []
    nuc_num = 0
    for lines in fasta_fh:
        if lines.startswith(">"):
            continue
        lines = lines.rstrip("\n")
        for i,letter in enumerate(lines):
            nuc_num = nuc_num + 1
            for juncts in junctions_list:
                frag1 = juncts[2]
                frag2 = juncts[3]
                lig_count_2 = juncts[4]

                if nuc_num == frag1:
                    print (nuc_num, letter, nuc_num-1, nuc_num+1, frag2, lig_count_2, sep ="\t", file = lig_connect_fh)
                    print (nuc_num, frag2, 1, lig_count_2, sep ="\t", file = lig_helix_fh)

                    #Writing a filtered connect file. Skip if first and last residue.
                    #Write to file if frag 1 nucleotide is G (3' of fragment 1 in ligation file)
                if (nuc_num == frag1 == 1) or (nuc_num == frag1 == max_num_num):
                        continue

                elif nuc_num == frag1 and letter == 'G':
                    print (nuc_num, letter, nuc_num-1, nuc_num+1, frag2, lig_count_2, sep ="\t", file = lig_connect_filtered_fh)
                    print (frag1, frag2, 1, lig_count_2, sep ="\t", file = lig_helix_filtered_fh)


    print ("\nWrote ligation junctions in a connect file and filtered connect file!\n")
    print ("\nWrote ligation junctions in a helix file and filtered helix file!\n")

    fasta_fh.close()
    lig_helix_fh.close()
    lig_connect_fh.close()
    lig_helix_filtered_fh.close()
    lig_connect_filtered_fh.close()


def junction_frequency(connect_file, junction_nuc, out_file):
    connect_fh = open(connect_file,'r')
    out_junc_nuc_fh = open(out_file+'_juncnuc_basepair.txt','w')

    print ("Position (connect)", "Nucleotide", "Base-paired Nucleotide", "Position (ligation)", "Frequency in Junctions",
                sep = "\t", file = out_junc_nuc_fh)

    #Identify base-paired nucleotides
    connect_fh.readline()
    for i,base in enumerate(connect_fh):
        base = base.rstrip("\n").split(" ")

        base_num = int(base[0])
        base_res = base[1]
        base_pair = base[4]

        #Get Frequency of each nucleotide in ligation
        for key, value in Counter(junction_nuc).items():
            if key == base_num:
                #Write a file with
                print (base_num, base_res, base_pair, key, value, sep ="\t", file = out_junc_nuc_fh)

    out_junc_nuc_fh.close()
    print ("\nWrote ligation vs Base pairs file!\n")


if __name__=="__main__":
    main()


"""
Compares ligation Junction to 3D distances obtained from PDB file. Also, adds base pairing info if available from .ct file
Ligation junction  from 'itx2scatterplot_strand.py script'. Sample line from file-
RNA5-8S5       	RNA5-8S5       	74     	108    	6
RNA5-8S5       	RNA5-8S5       	136    	56     	2

3D distances from Pymol_distance.py script. Sample line from file-
resi1(4ug0_5.8S)       	resi2(4ug0_5.8S)       	pdb_id_1       	pdb_id_2       	distance
2      	2      	82762  	82762  	0.0
2      	3      	82762  	82785  	5.787150383

Base pairing info from connect file-
73 ENERGY =     -17.50    S.cerevisiae_tRNA-PHE
 1 G       0    2   72    1
 2 C       1    3   71    2
 3 G       2    4   70    3
 4 G       3    5   69    4
 The connect format is column based. The first column specified the sequence index, starting at one.
 Columns 3, 4, and 6 redundantly give sequence indices (plus/minus one).
 The second column contains the base in one-letter notation.
 Column 4 specifies the pairing partner of this base if it involved in a base pair.
 If the base is unpaired, this column is zero.
"""
