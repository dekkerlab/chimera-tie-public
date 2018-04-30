#!/bin/env python3
import glob
import argparse
import gzip
import os

def get_arguments():
    parser = argparse.ArgumentParser(description="""
    Given the matrices m1 and m2,
    it sets all the pixels (i,j) in m1 to 0 if
    the pixel (i,j) in m2 > 0.

    TODO
    Write the naming conventions
    as this script heavily relies on file naming conventions
    with underscores
    EX:
    Rep1_minus_XIST_10nt_D.matrix.gz
    Rep1_plus_XIST_10nt_D.matrix.gz
    """)
    parser.add_argument("-i" ,
                        help = """Input Folder""" ,
                        required = False ,
                        type = str)
    parser.add_argument("-o" ,
                        help = """Output Folder""" ,
                        required = False ,
                        type = str)

    arguments = parser.parse_args()
    return arguments

def get_file_pairs(matrix_folder):
    plus_files = sorted(list(glob.glob(matrix_folder + "/*_plus_*matrix.gz") ))
    minus_files = sorted(list(glob.glob(matrix_folder + "/*_minus_*matrix.gz") ))

    if len(plus_files) != len(minus_files):
        exit("The number of plus files are different from the minus files")

    return zip(plus_files, minus_files)

def filter_matrix( file_pair, output_folder ):
    output_file = "?????"
    transcript_name = os.path.basename(file_pair[0]).split("_")[2]
    output_file = os.path.join(output_folder,
                               transcript_name+".filtered.matrix.gz")


    with gzip.open(file_pair[0], "rt") as plus_stream,\
         gzip.open(file_pair[1], "rt") as minus_stream,\
         gzip.open(output_file, "wt") as output_stream:

         plus_first_line = plus_stream.readline().strip()
         minus_first_line = minus_stream.readline().strip()

         while plus_first_line[0] == "#":
              plus_first_line = plus_stream.readline().strip()
         while minus_first_line[0] == "#":
              minus_first_line = minus_stream.readline().strip()

         plus_dims = plus_first_line.strip().split()[0].split("x")
         minus_dims = minus_first_line.strip().split()[0].split("x")

         print(plus_first_line, file=output_stream)

         if (plus_dims[0] != minus_dims[0]) or\
            (plus_dims[1] != minus_dims[1]):
            exit("Mismatch in dimensions of" +  str(file_pair) )

         # Numpy conflicts due to python version in the current cluster setup
         # so we are doing the filtering in a slow way
         # by going through the elements one by one
         for x in plus_stream:
             plus_line = x.strip().split()
             minus_line = minus_stream.readline().strip().split()
             output_array = list()
             for index, interaction in enumerate(plus_line[1:]):
                 if float(minus_line[index + 1]) > 0:
                     output_array.append(0)
                 else:
                     output_array.append(interaction)

             this_line = plus_line[0] + "\t" + "\t".join(map(str, output_array))
             print(this_line, file = output_stream)


def main():
    arguments = get_arguments()
    file_pairs = get_file_pairs( arguments.i )
    if not os.path.isdir(arguments.o):
        os.makedirs(arguments.o)

    for k in file_pairs:
        filter_matrix(k, arguments.o)


if __name__ == "__main__":
    main()
