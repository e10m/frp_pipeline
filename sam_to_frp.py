#!/home/ethanmach1998/anaconda3/envs/snakemake/bin/python

import os
import sys
import re
import matplotlib.pyplot as plt
import statistics


def parse_sam(sam_file):
    """
    This Python script accepts SAM files and parses it, returning a tuple containing two lists and a variable

    :param sam_file: Input SAM alignment file that includes headers and MD lines
    :return: tuple ([pos_list], [percent_identities], taxon_id)
    """
    # initialize variables / data structures
    pos_list = []
    perc_identity_list = []

    with open(sam_file, 'r') as file:
        for line in file:
            # obtain ref_genome name
            if line.startswith("@SQ"):
                ref_genome = re.search(r'SN:(\S+)', line).group(1)
                continue

            # skip header lines
            elif line.startswith("@"):
                continue

            # calculate total number of matches from CIGAR line
            total_cigar_matches = 0
            actual_matches = 0
            deletions = 0

            cigar_line = re.search(r'\*|([0-9]+[MIDNSHPX=])+', line).group(0)
            match_list = re.findall(r'(\d+)M', cigar_line)

            for match in match_list:
                match = int(match)
                total_cigar_matches += match

            # calculate stats from MD line
            md_line = re.search(r'(MD):(Z):(([0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*))', line).group(3)
            md_matches = re.findall(r'(\d+)', md_line)
            md_deletions = re.findall(r'\^([A-Z]+)', md_line)

            # calculate actual # of matches
            for match in md_matches:
                match = int(match)
                actual_matches += match

            # calculate # of deletions
            for deletion in md_deletions:
                deletions += len(deletion)

            # calculate percent identity and append it to a list
            percent_identity = 100 * (actual_matches/(total_cigar_matches + deletions))
            perc_identity_list.append(percent_identity)

            # search for and append pos starting position
            pos = re.search(r'(\d+)\t(\d+)\t(\*|([0-9]+[MIDNSHPX=])+)', line)
            pos_list.append(pos.group(1))

    # returns a tuple of ([positions], [percent identities])
    return pos_list, perc_identity_list, ref_genome


def generate_frp(data, file):
    """
    This function takes parsed SAM data and the file path name, generating a fragment recruitment plot (FRP) via
    matplotlib

    :param data: Parsed sam data in a tuple ([pos_list], [percent_identities], taxon_id)
    :param file: file path name (str)
    :return: Fragment Recruitment Plot image (.png)
    """
    # store data
    positions = data[0]
    percent_identities = data[1]
    ref_genome = data[2]

    # parse file name
    file = re.search(r'(/\w+/\w+/\w+/\w+)/(\d+)/\d+_([a-z]+\d*)', file)

    taxon_id = file.group(2)
    alignment_method = file.group(3)
    hmmc_dir = file.group(1)

    print("The SAM data has been loaded in")

    # generate plot
    plt.scatter(positions, percent_identities, color='blue', marker='o', s=1, alpha=0.5)
    print("The fragment recruitment plot has been generated")

    # label plot
    plt.title(f'FRP Against Reference Genome {ref_genome} Using {alignment_method}', fontsize=10)
    plt.xlabel('Position on Reference Genome')
    plt.ylabel('Percent Identity')
    plt.ylim(80, 100)
    print("The chart has been labeled")


    # write/output plot as .png image
    plt.savefig(f'{hmmc_dir}/{taxon_id}/{taxon_id}_{alignment_method}_frp.png', format='png', dpi=300)
    print("The plot image file has been created")


if __name__ == '__main__':
    file_path = sys.argv[1]

    # open and parse the input SAM file
    with open(file_path, "r") as input_file:
        parsed_sam_data = parse_sam(file_path)

    print("The SAM file has been parsed successfully")

    # plot a fragment recruitment plot using the parsed sam data
    generate_frp(parsed_sam_data, file_path)