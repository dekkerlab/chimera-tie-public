#!/bin/env python3
import glob
import argparse
import gzip
import os

def get_arguments():
    parser = argparse.ArgumentParser(description="""
    Given a folder fo *matrix.gz files,
    compute the total interactions and normalize it by matrix size.
    The output is a tab separated list of values as follows

    #Transcript_Name   Total_Interactions Dimension Total_Interactions/(Dimension/2)*(Dimension - 1)
    """)
    parser.add_argument("-i" ,
                        help = """Input Folder""" ,
                        required = False ,
                        type = str)
    parser.add_argument("-o" ,
                        help = """Output File""" ,
                        required = False ,
                        type = str)

    arguments = parser.parse_args()
    return arguments

# Note that we are ignoring the first diagonal
# and summing the values in the upper triangle
def get_interactions(matrix_file):
    dimension = -1
    total_interactions = 0
    row_counter = 1

    transcript = os.path.basename(matrix_file).split(".matrix.gz")[0]

    with gzip.open(matrix_file, "rt") as input_stream:
        first_line = input_stream.readline().strip()

        while first_line[0] == "#":
             first_line = input_stream.readline().strip()

        dimension = int(first_line.split()[0].split("x")[0])

        for this_line in input_stream:
            contents = this_line.strip().split()
            upper_diagonal = contents[row_counter + 1 :]
            this_sum = sum(map(float, upper_diagonal))
            total_interactions += this_sum
            row_counter += 1
        normalized_sum = total_interactions / ( (dimension/2)*(dimension -1)  )

    return [ transcript, total_interactions, dimension, normalized_sum ]


def get_files(matrix_folder):
    return list(glob.glob(matrix_folder + "/*matrix.gz"))



def main():
    arguments = get_arguments()
    matrix_files = get_files(arguments.i)
    output_array = list()

    for f in matrix_files:
        interactions = get_interactions(f)
        output_array.append(interactions)

        with open(arguments.o, "w") as output_stream:
            print("\t".join(["#Transcript_Name",
            "Total_Interactions", "Dimension", "Normalized_Interactions"]),
            file = output_stream)
            for e in output_array:
                print("\t".join( map(str, e) ), file = output_stream)


if __name__ == "__main__":
        main()
