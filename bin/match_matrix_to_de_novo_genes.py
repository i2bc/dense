#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes as input a match_matrix file and a TRG_before_filtering file, both in tsv format. 

The match_matrix file should contain information about the matches, while the TRG_before_filtering file should contain information about the parent genes and their associated CDS (Coding DNA Sequences).

The script supports three strategies for defining de novo genes, which can be specified using the --strategy argument. The strategies are as follows:

1. If synteny is required, the pattern is "gS*". If not, the pattern is either "gS*" or "gNS*".
2. If synteny is required, the pattern is either "gS*" or "gS". If not, the pattern can be "gS*", "gS", "gNS*", or "gNS".
3. If the second column is "CDS" and "CDS" is not in any other column, the pattern is "gS*" if synteny is required. If not, the pattern is either "gS*" or "gNS*".

The script outputs a tsv file containing the de novo genes. The path to the output file can be specified using the -o or --output argument.

Usage:
    python match_matrix_to_de_novo_genes.py match_matrix.tsv TRG_before_filtering.tsv --strategy 1 --synteny True -o output.tsv
"""
import argparse
import csv
import sys

# Create the parser
parser = argparse.ArgumentParser(description='Match matrix to de novo genes')

# Add the arguments
parser.add_argument('match_matrix', type=str, help='The path to the match_matrix tsv file')
parser.add_argument('TRG_before_filtering', type=str, help='The path to the TRG_before_filtering tsv file')
parser.add_argument('--strategy', type=int, choices=[1, 2, 3], default=1, help='The strategy to define de novo genes')
parser.add_argument('--synteny', type=str, help='Whether or not to require that non-coding matches are in synteny')
parser.add_argument('--num_outgroups', type=int, help='The required number of outgroups with a non-coding match for the strategy 1')
parser.add_argument('-o', '--output', type=str, help='The path to the output tsv file')

# Parse the arguments
args = parser.parse_args()


#### FILTER THE MATCH MATRIX ####
de_novo_CDS = []
not_de_novo_CDS = []

# The maximum number of outgroups is the number of columns in the match matrix minus 2 (first and second columns)
with open(args.match_matrix, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    num_columns = len(header)
    max_num_outgroups = num_columns - 2
# Check if the number of outgroups is valid
if args.num_outgroups > max_num_outgroups or args.num_outgroups <= 0:
    sys.exit(f"Invalid number of outgroups. The maximum number of outgroups is {max_num_outgroups}")

# Open the match_matrix file
with open(args.match_matrix, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # Skip the header

    # Process each line
    for line in reader:
        # Remove the eventual terminal "_elongated" from the CDS name
        if line[0].endswith("_elongated"):
            line[0] = line[0].replace("_elongated", "")

        isDeNovo = False
        # Possible values : "noM", "CDS", "gS", "gNS", "gNA". The last 3 values can be followed by a number if they are outgroups (e.i. after the last CDS match (strategy 1))
        if args.strategy == 1 : # Must have at least n outgroups with a non-coding match
            line_outgroup_synteny = set(column[2:] for column in line[2:] if column.startswith("gS") and column[2:].isdigit())

            if args.synteny == "true":
                # If their is n different values of pattern "gS[0-9]" in the line (starting from the third column)
                if len(line_outgroup_synteny) >= args.num_outgroups:
                    isDeNovo = True
            else:
                # starts.with("gN") can be either "gNS" or "gNA"
                line_outgroup_non_synteny = set(column[3:] for column in line[2:] if column.startswith("gN") and column[3:].isdigit())
                if len(line_outgroup_synteny) + len(line_outgroup_non_synteny) >= args.num_outgroups:
                    isDeNovo = True

        elif args.strategy == 2 :
            if args.synteny == "true":
                if any(column.startswith("gS") for column in line[2:]):
                    isDeNovo = True
            else:
                if any(column.startswith("gS") or column.startswith("gN") for column in line[2:]):
                    isDeNovo = True

        elif args.strategy == 3 :
            if line[1] == "CDS" and not any("CDS" in column for column in line[2:]): # Check that their is one single CDS match (query)
                # Then same as strategy 2
                if args.synteny == "true":
                    if any(column.startswith("gS") for column in line[2:]):
                        isDeNovo = True
                else:
                    if any(column.startswith("gS") or column.startswith("gN") for column in line[2:]):
                        isDeNovo = True

        if isDeNovo:
            # Store the CDS as de novo
            de_novo_CDS.append(line[0])
        else :
            # Store the CDS as not de novo
            not_de_novo_CDS.append(line[0])


# Write a file named "TRG_after_strategy_filtering.txt"
with open("TRG_after_strategy_filtering.txt", 'w') as f:
    for CDS in de_novo_CDS:
        f.write(CDS + "\n")


#### ISOFORMS CONTROL ####
# Read the TRG_before_filtering file
with open(args.TRG_before_filtering, 'r') as f:
    reader = csv.reader(f, delimiter='\t')

    # Create a dictionary of parent genes and their CDS
    parent_genes = {}
    for line in reader:
        CDS, parent_gene = line
        if parent_gene in parent_genes:
            parent_genes[parent_gene].append(CDS)
        else:
            parent_genes[parent_gene] = [CDS]

# Check which parent genes have all their CDS in the de_novo_CDS list
de_novo_genes = {}
discarded_de_novo_genes = {}
for parent_gene, CDS_list in parent_genes.items():
    if any(CDS in de_novo_CDS for CDS in CDS_list) and any(CDS in not_de_novo_CDS for CDS in CDS_list) :
        discarded_de_novo_genes[parent_gene] = CDS_list
    elif all(CDS in de_novo_CDS for CDS in CDS_list) :
        de_novo_genes[parent_gene] = CDS_list

# If discarded_de_novo_genes is not empty, write a file named "TRG_discarded_because_of_isoforms.tsv"
if discarded_de_novo_genes:
    with open("TRG_discarded_because_of_isoforms.tsv", 'w') as f:
        f.write("gene" + "\t" + "CDS" + "\n")
        for parent_gene, CDS_list in discarded_de_novo_genes.items():
            for CDS in CDS_list:
                f.write(parent_gene + "\t" + CDS + "\n")

#### OUTPUT ####
# Write the output file
with open(args.output, 'w') as f:
    f.write("gene" + "\t" + "CDS" + "\n")
    for parent_gene, CDS_list in de_novo_genes.items():
        for CDS in CDS_list:
            f.write(parent_gene + "\t" + CDS + "\n")