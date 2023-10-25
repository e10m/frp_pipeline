#!/usr/bin/python3

import os
import sys
import re


def combine_contigs(fasta_file):
    """
    This Python script combines fragmented reference genomes via an artificially generated linker of nucleotides (XXXXX)
    :param file: FASTA file expected (.fna, .faa, .fasta, etc.)
    :return: A new file with the contigs linked
    """
    # initialize variables
    contig_list = []
    header_line = None
    new_fasta_list = []

    with open(fasta_file, 'r') as file:
        for line in file:
            # save header for writing out to file later
            if line.startswith(">") and header_line is None:
                # regex search
                header_line = re.search(r'(>\w+\d+\.*\d*)\s(\w+\s\w+)', line)
                header_line = header_line.group(0)
            # add linker chain of nucleotides
            elif line.startswith(">") and header_line is not None:
                contig_list.append("XXXXX")
            # save actual sequence
            else:
                line = line.rsplit()
                contig_list.append(line[0])

    # combine the sequences together
    combined_contig = "".join(contig_list)

    new_fasta_list.append(header_line + ' linked genome\n')

    # format combined sequence into 80 characters per line
    for i in range(0, len(combined_contig), 80):
        new_line = combined_contig[i:i+80]
        new_fasta_list.append(new_line+'\n')

    # combine sequences together again
    new_fasta = "".join(new_fasta_list)

    return new_fasta


if __name__ == '__main__':
    file_path = sys.argv[1]

    # error handling for if file path exists
    if not os.path.exists(file_path):
        print("No file is found at {}".format(file_path))
        sys.exit(1)

    # open and manipulate file
    with open(file_path, "r") as input_file:
        output_fasta = combine_contigs(file_path)

    output_path = file_path.rsplit('.', 1)[0] + "_combined.fna"

    # write to new file
    with open(output_path, 'w') as output_file:
        output_file.write(output_fasta)

    print("The new file has been written to: {}".format(output_path))
