#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:19:21 2022

@author: paul.roginski

The idea behind this script is that all possible most recent common ancestors 
of : a given species on the first hand, and its neighbors on the other hand,
are simply all the (node-represented) ancestors of the given species.
AND, all the ancestors-nodes of a given species, can be order by counting the 
number of edges between them and the root.
"""


import argparse
import dendropy
import operator
import csv
import sys


# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-focal", required=True, help="name of the focal species")
parser.add_argument("-out", required=False, help="output file name")
args = parser.parse_args()

  
# newick file for the phylogeny tree
tree_file = args.tree
tree = dendropy.Tree.get(path=tree_file,
                         schema='newick',
                         preserve_underscores=True,
                         rooting='force-rooted')



#print(tree.as_ascii_plot())
    
   
# Retrieve in the taxon_namespace attribute of the tree object, the element 
# corresponding with the focal species.
focal_name = args.focal
print("focal name : {}".format(focal_name))
taxon_labels = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]
print("names : {}".format(taxon_labels))




# TIME OF DIVERGENCE
# Build a distance matrix object
pdc = tree.phylogenetic_distance_matrix()


# Get the focal taxon
index_list = [ i for i in range(0,len(tree.taxon_namespace)) if tree.taxon_namespace[i].label == focal_name ]
if len(index_list) > 0 :
	focal = tree.taxon_namespace[index_list[0]]
else : print("'{}' could not be found in the provided tree.".format(focal_name)) ; sys.exit()


# Build a dictionnary with species as keys and distance to the focal species as values
# /2 to turn patristic distance into time of divergence
distance_to_focal = { t.label: pdc(focal, t)/2 for t in tree.taxon_namespace}
# Sorted version
distance_to_focal = dict(sorted(distance_to_focal.items(), key=operator.itemgetter(1)))




# ROOT DISTANCE
# Set every edge length to 1. By doing so, we are sure the calc_node_root_distances can be used
# Besides, we do not care about the actual value of each edge length in this case.
for node in tree.preorder_node_iter():
    node.edge.length = 1

# Adds attribute “root_distance” to each node, with value set to the sum of 
# edge lengths from the node to the root. Returns list of distances. 
tree.calc_node_root_distances(return_leaf_distances_only=False) 

# For each taxon provided in names, get its most recent common ancestor with 
# the focal species (= a node), and get its distance to the tree root.
root_distance = { name:tree.mrca(taxon_labels=[focal_name, name]).root_distance for name in taxon_labels if name != focal_name }
root_distance[focal_name]=max(root_distance.values()) +1
print("root_distance : {}".format(root_distance))


final_dic = [ {"Taxon" : key, "Time of divergence (MYA)" : distance_to_focal[key], "Root distance" : root_distance[key]} for key in distance_to_focal ]

fields = ["Taxon", "Time of divergence (MYA)", "Root distance"]

# Open the CSV file with write permission
with open(args.out, "w", newline="") as csvfile:

    writer = csv.DictWriter(csvfile, fieldnames=fields, delimiter = "\t", lineterminator="\n")
    writer.writeheader()
    for row in final_dic:
        writer.writerow(row)
