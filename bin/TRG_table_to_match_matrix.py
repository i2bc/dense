#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:35:34 2023

@author: Proginski
"""

# This script takes a TRG_table and returns a match matrix where genomes are columns and CDS (TRG) are rows
# Columns (genomes) are sorted by distance to the focal genome if a tree is provided.

import argparse
from collections import defaultdict
import dendropy
import sys
import os

# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("input", help="TRG_table")
parser.add_argument("tree", help="newick tree")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


#### MATCHES COLLECTION ####
# Build a dictionary with CDS as keys and as values as sub-dictionary with :
# CDS : the list of neighbors having a CDS match
# genome_syntNA : the list of neighbors having a genome (synteny not checked)
# genome_synt : the list of neighbors having a genome match in synteny
# genome_nosynt : the list of neighbors having a genome match not in synteny
match_dic = defaultdict(lambda: {'CDS':[], 'genome_syntNA':[], 'genome_synt':[], 'genome_nosynt':[]})
neighbors = []

# Open the file once to get the focal genome
with open(args.input) as input_file:
    # Skip the first line
    next(input_file)

    # Get focal genome from the second line
    second_line = next(input_file).strip().split('\t')
    focal_genome = second_line[0]
    neighbors.append(focal_genome)

# Open the file a second time to process the lines
with open(args.input) as input_file:
    # Skip the first line
    next(input_file)

    # Process the rest of the lines
    for line in input_file:
        line = line.strip().split('\t')
        cds = line[1]
        neighbor = line[2]
        match_type = line[4]

        # Record unique neighbors
        if neighbor not in neighbors:
            neighbors.append(neighbor)

        if match_type == "CDS":
            match_dic[cds]['CDS'].append(neighbor)
        elif match_type == "genome":
            isSynt = line[7]
            if isSynt == "NA":
                match_dic[cds]['genome_syntNA'].append(neighbor)
            elif isSynt == "True":
                match_dic[cds]['genome_synt'].append(neighbor)
            elif isSynt == "False":
                match_dic[cds]['genome_nosynt'].append(neighbor)


#### TREE INFORMATION ####
# Open the tree file to order the neighbors with respect to the focal genome.
# To do so, get the distance between the focal genome and each neighbor.
# Two of more neighbors of the same "outgroup" should have the same distance to the focal genome.
# The classical 'time of divergence' may lead to small differences if branches are not of equal length.
# For example, Pan paniscus and Pan troglodytes form a possible "outgroup" in the human-and-its-neighbors tree,
# but their calculated time of divergence with human may not be the same.
# So just get the rank of the most recent common ancestor (MRCA) node.

# Load the tree (if non-empty)
if os.stat(args.tree).st_size > 0:
    tree = dendropy.Tree.get(path=args.tree,
                            schema='newick',
                            preserve_underscores=True,
                            rooting='force-rooted')
    
    # Get the list of taxon labels (e.g. genome names)
    taxon_labels = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]

    # Build a distance matrix object
    pdc = tree.phylogenetic_distance_matrix()

    # Get the focal taxon object
    focal_index = [ i for i in range(0,len(tree.taxon_namespace)) if tree.taxon_namespace[i].label == focal_genome ]
    try :
        focal = tree.taxon_namespace[focal_index[0]]
    except :
        sys.exit(f"{focal_genome} could not be found in the provided tree.")

    # Set every edge length to 1. By doing so, we are sure the calc_node_root_distances can be used
    # Besides, we do not care about the actual value of each edge length in this case.
    for node in tree.preorder_node_iter():
        node.edge.length = 1

    # Adds attribute “root_distance” to each node, with value set to the sum of 
    # edge lengths from the node to the root. Returns list of distances. 
    tree.calc_node_root_distances(return_leaf_distances_only=False) 

    # For each genome in taxon_labels, get its most recent common ancestor with 
    # the focal genome (= a node), and get its distance to the tree root.
    root_distance = { name:int(tree.mrca(taxon_labels=[focal_genome, name]).root_distance) for name in taxon_labels if name != focal_genome }
    root_distance[focal_genome]=max(root_distance.values()) +1
    # Reverse the dictionary to make further comparisons more intuitive (focal is 0, its closest neighbor is 1, etc.)
    max_distance = max(root_distance.values())
    focal_distance = { genome: max_distance - root_distance[genome] for genome in root_distance }

    # List of genomes sorted by distance to the focal genome
    sorted_genomes = sorted(focal_distance, key=focal_distance.get)


    #### MATCH MATRIX ####
    with open(args.output, 'w') as output_file :
        # Write the header
        output_file.write("CDS\t{}\n".format('\t'.join(sorted_genomes)))
        # Write the rest of the matrix, for each CDS
        for cds in match_dic :
            line = [cds]
            # Get the maximum distance for a CDS match
            # Genome matches that are farther than this distance are "outgroup" matches
            max_CDS_distance = max(focal_distance[genome] for genome in sorted_genomes if genome in match_dic[cds]['CDS'])
            # Write the line
            for genome in sorted_genomes :
                if genome in match_dic[cds]['CDS']:
                    line.append("CDS")
                else : # 'genome' match
                    suffix=" "
                    if focal_distance[genome] > max_CDS_distance :
                        # For matches in an outgroup species, the suffix becomes a digit indicating the rank of the outgroup
                        suffix=str(focal_distance[genome]-max_CDS_distance)
                    if genome in match_dic[cds]['genome_syntNA']:
                        line.append("gNA"+suffix)
                    elif genome in match_dic[cds]['genome_synt']:
                        line.append("gS"+suffix)
                    elif genome in match_dic[cds]['genome_nosynt']:
                        line.append("gNS"+suffix)
                    else :
                        line.append("noM")
            output_file.write("{}\n".format('\t'.join(line)))
            
else : # No tree provided
    #### MATCH MATRIX ####
    with open(args.output, 'w') as output_file :
        # Write the header
        output_file.write("CDS\t{}\n".format('\t'.join(neighbors)))
        # Write the rest of the matrix, for each CDS
        for cds in match_dic :
            line = [cds]
            for genome in neighbors :
                if genome in match_dic[cds]['CDS']:
                    line.append("CDS")
                else : # 'genome' match
                    suffix=" "
                    if genome in match_dic[cds]['genome_syntNA']:
                        line.append("gNA"+suffix)
                    elif genome in match_dic[cds]['genome_synt']:
                        line.append("gS"+suffix)
                    elif genome in match_dic[cds]['genome_nosynt']:
                        line.append("gNS"+suffix)
                    else :
                        line.append("noM")
            output_file.write("{}\n".format('\t'.join(line)))
