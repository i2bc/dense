#!/usr/bin/env python3

import sys
import os
import subprocess
import os
import re

def trg_node_from_trg_rank(taxdump, taxid, TRG_rank):
    """
    Retrieve the TRG_node threshold from the given TRG_rank.

    Args:
        taxdump (str): The path to the taxdump directory.
        taxid (str): The taxonomic ID.
        TRG_rank (str): The target rank.

    Returns:
        None

    Raises:
        SystemExit: If the TRG_rank could not be found in the phylogeny.
    """
    print(f"Retrieving TRG_node from the TRG_rank '{TRG_rank}'.")

    original_taxid = taxid
    rank = ""
    counter = 0

    while rank != TRG_rank and counter < 100:
        old_taxid = taxid
        with open(f"{taxdump}/nodes.dmp") as file:
            for line in file:
                if line.startswith(f"{taxid}\t"):
                    taxid = line.strip().split("\t|\t")[1]
                    rank = line.strip().split("\t|\t")[2]
                    print(f"taxid: {taxid}, rank: {rank}")

                    break
        counter += 1

    print("")

    if rank != "no rank" and counter != 100:
        with open(f"{taxdump}/fullnamelineage.dmp") as file:
            for line in file:
                if line.startswith(f"{old_taxid}\t"):
                    TRG_node = line.strip().split("\t|\t")[1]
    else:
        print(f"{TRG_rank} could not be found in the phylogeny of taxid {original_taxid}.")
        sys.exit(1)

    return TRG_node

def get_trg_nodes(taxdump, taxid, TRG_node, TRG_rank):
    """
    Filter the genera output based on the provided parameters.

    Args:
        taxdump (str): Path to the taxdump directory.
        taxid (str): The taxid of the focal genome to filter the genera output.
        TRG_node (str): The target node in the taxonomic hierarchy. If not provided, it will be determined based on the taxid, and the TRG_rank.
        TRG_rank (str): The rank of the target node in the taxonomic hierarchy.

    Returns:
        list: List of TRG nodes that will be used to define TRGs.
    """
    if taxid == "EMPTY":
        print("The focal taxid could not be found!")
        sys.exit(1)

    lineage_file = os.path.join(taxdump, 'fullnamelineage.dmp')
    if not( os.path.isfile(lineage_file) and os.path.getsize(lineage_file) ) > 0:
        print(f"The file {lineage_file} does not exist or is empty.")
        sys.exit(1)

    # If TRG_node is not provided, use TRG_rank
    if TRG_node is None:
        TRG_node = trg_node_from_trg_rank(taxdump, taxid, TRG_rank)
    else:
        TRG_node = TRG_node.strip()

    print(f"TRG_node is now set to '{TRG_node}'.")

    # Identify the TRG nodes
    print(f"Looking for taxid {taxid} in {lineage_file}...")
    # Continue with the rest of the code here
    TRG_nodes = []
    with open(lineage_file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(f"{taxid}\t"):
                print(f"Found :\n{line.strip()}")
                taxid_node = line.strip().split("\t|\t")[1]
                if taxid_node == TRG_node:
                    TRG_nodes = [taxid_node]
                else :
                    lineage = line.strip().split("\t|\t")[2]
                    nodes = lineage.strip("; \t|").split("; ")
                    found_TRG_node = False
                    for node in nodes:
                        if node == TRG_node:
                            found_TRG_node = True
                        if found_TRG_node:
                            TRG_nodes.append(node)
                    if found_TRG_node:
                        TRG_nodes.append(taxid_node)
                break

    if not TRG_nodes:
        raise ValueError(f"{TRG_node} could not be found. Cannot proceed with filtering.")
    
    # Print TRG_nodes
    print("\nNodes that will be used to define TRGs:")
    for node in TRG_nodes:
        print(node)

    return TRG_nodes

def filter_genera_output(taxdump, taxid, TRG_node, TRG_rank, focal_mRNA_to_gene, genera_out):
    """
    Filter the genEra output based on the provided parameters.

    The genEra output is a tab delimited file with the following header :
    '#gene	phylostratum	rank	taxonomic_representativeness'

    The first column contains mRNA/CDS/isoform/protein names and the second column, a phylostratum (e.i. a 'node').

    The function :
    - defines the 'oldest' node to consider TRGs (e.i. the 'TRG_node'),
    - retrieves all the nodes from the TRG_node to present (e.i 'TRG_nodes'),
    - adds a terminal column to the genEra output with the gene names
    - write a file 'TRG_CDS.txt' with the mRNA of genes whose all mRNA are young.
    """
    TRG_nodes = get_trg_nodes(taxdump, taxid, TRG_node, TRG_rank)

    # From the focal_mRNA_to_gene mapping file, build a dictionary with mRNA as key and mRNA as gene.
    # Also do the reverse.
    mapping_genes = {}
    with open(focal_mRNA_to_gene) as f:
        for line in f:
            mRNA, gene = line.strip().split("\t")
            mapping_genes[mRNA] = gene
    mapping_mRNAs = {}
    with open(focal_mRNA_to_gene) as f:
        for line in f:
            mRNA, gene = line.strip().split("\t")
            if gene in mapping_mRNAs:
                mapping_mRNAs[gene].append(mRNA)
            else:
                mapping_mRNAs[gene] = [mRNA]

    # Write a new version of the genEra output with an additional column with the gene name at the end.
    genera_mRNAs = set()
    with open(genera_out) as f, open('genEra_output_with_genes.tsv', 'w') as out_file, open('missing_in_mapping.txt', 'w') as missing_file:
        for line in f:
            cols = line.strip().split("\t")
            mRNA = cols[0]
            if mRNA not in genera_mRNAs:
                genera_mRNAs.add(mRNA)
            if mRNA in mapping_genes:
                out_file.write(f"{line.strip()}\t{mapping_genes[mRNA]}\n")
            elif mRNA != "#gene":
                print(f"WARNING: {mRNA}'s gene could not be found in the mapping file.", file=sys.stderr)
                missing_file.write(f"{mRNA}\n")

    # Write a file with the mRNA that are in the mapping file but not in the genEra output.
    with open('missing_in_genEra_output.txt', 'w') as f:
        for gene in mapping_mRNAs:
            for mRNA in mapping_mRNAs[gene]:
                if mRNA not in genera_mRNAs:
                    print(f"WARNING: {mRNA} is not in the genEra output.", file=sys.stderr)
                    f.write(f"{mRNA}\n")

    # Udpate mapping_mRNAs to only keep the mRNAs that are in the genEra output.
    for gene in mapping_mRNAs:
        mapping_mRNAs[gene] = list(set(mapping_mRNAs[gene]).intersection(genera_mRNAs))

    # Build a set of genes that have at least one isoform whose 'age' is in the TRG_nodes,
    # and a set of genes that have at least one  isoform older than the TRG_nodes.
    genes_with_young = set()
    genes_with_old = set()
    with open('genEra_output_with_genes.tsv') as f:
        for line in f:
            cols = line.strip().split("\t")
            gene = cols[-1]
            node = cols[1]
            if node in TRG_nodes :
                if gene not in genes_with_young:
                    genes_with_young.add(gene)
            else :
                if gene not in genes_with_old:
                    genes_with_old.add(gene)

    # TRGs are the genes that only have young isoforms.
    TRGs = list(genes_with_young.difference(genes_with_old))

    # Main output file
    with open('TRG_CDS.txt', 'w') as f:
        for TRG in TRGs :
            for mRNA in mapping_mRNAs[TRG] :
                f.write(f"{mRNA}\n")

    # For information
    with open('TRG_nodes.txt', 'w') as f:
        for node in TRG_nodes :
            f.write(f"{node}\n")
    with open('TRG_before_isoform_correction.txt', 'w') as f:
        for gene in genes_with_young :
            f.write(f"{gene}\n")

if __name__ == "__main__":
    taxdump = sys.argv[1]
    taxid = sys.argv[2]
    TRG_node = sys.argv[3] if sys.argv[3] != "null" else None
    TRG_rank = sys.argv[4]
    focal_mRNA_to_gene = sys.argv[5]
    genera_out = sys.argv[6]

    filter_genera_output(taxdump, taxid, TRG_node, TRG_rank, focal_mRNA_to_gene, genera_out)
