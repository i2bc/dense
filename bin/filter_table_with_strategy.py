#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:55:42 2023

@author: paul.roginski
"""

import argparse
import re

def filter_table(table_path, strategy, out_path, check_synteny=False):
    """
    Filters a table of TRG-hit relationships based on a specified strategy.
    
    Args:
        table_path (str): The path to the table of TRG-hit relationships.
        strategy (int): The strategy to apply for filtering.
            1: Looks for a non-coding hit in an outgroup species.
            2: Looks for a non-coding hit.
            3: Looks for orphan genes.
        out_path (str): The path to the output file.
    
    Raises:
        ValueError: If the specified strategy is not between 1 and 3.
    """

    if not strategy in range(1,4):
        raise ValueError("'--strategy' can only be 1, 2 or 3.")
        
    if strategy == 1:
        table_dic={}
        # Store the root distance of CDS and genome hits for each TRG.
        with open(table_path, 'r') as table:
            next(table)
            for line in table:
                line=line.strip().split('\t')
                
                trg=line[1] # name of the TRG
                hit_type=line[4] # CDS or genome
                root_distance=float(line[6]) # the smallest, the farest from the focal
                
                if trg not in table_dic:
                    table_dic[trg]={
                            'CDS_hits_root_distance' : [],
                            'genome_hits_root_distance' : []
                            }
                if hit_type == "CDS":
                    table_dic[trg]['CDS_hits_root_distance'].append(root_distance)
                elif hit_type == "genome":
                    if check_synteny == True:
                        if line[7] == "True" :# Boolean
                            table_dic[trg]['genome_hits_root_distance'].append(root_distance)
                    else:
                        table_dic[trg]['genome_hits_root_distance'].append(root_distance)
                        
        with open(out_path, 'w') as out:
            for trg in table_dic:
                farest_CDS=float('inf')
                farest_genome=float('inf')
                
                if len(table_dic[trg]['CDS_hits_root_distance']) > 0 :
                    farest_CDS=min(table_dic[trg]['CDS_hits_root_distance'])
                if len(table_dic[trg]['genome_hits_root_distance']) > 0 :
                    farest_genome=min(table_dic[trg]['genome_hits_root_distance'])
                    
                # If the farest genome hits is farer than the farest CDS hit
                if farest_genome < farest_CDS:
                    trg=re.sub(r"_elongated.*", "", trg)
                    out.write(f"{trg}\n")
                    
    elif strategy == 2:
        table_dic={}
        # Store the root distance of CDS and genome hits for each TRG.
        with open(table_path, 'r') as table:
            next(table)
            for line in table:
                line=line.strip().split('\t')
                
                trg=line[1] # name of the TRG
                hit_type=line[4] # CDS or genome
                neighbor=line[2] # the smallest, the farest from the focal
                
                if trg not in table_dic:
                    table_dic[trg]={
                            'CDS_hits' : [],
                            'genome_hits' : []
                            }
                if hit_type == "CDS":
                    table_dic[trg]['CDS_hits'].append(neighbor)
                elif hit_type == "genome":
                    if check_synteny == True:
                        if line[7] == "True" :# Boolean
                            table_dic[trg]['genome_hits'].append(neighbor)
                    else:
                        table_dic[trg]['genome_hits'].append(neighbor)
                        
        with open(out_path, 'w') as out:
            for trg in table_dic:
                non_coding_hits = [ neighbor for neighbor in table_dic[trg]['genome_hits'] if neighbor not in table_dic[trg]['CDS_hits'] ]
                # If there is a non-coding hit
                if len(non_coding_hits) > 0 :
                    trg=re.sub(r"_elongated.*", "", trg)
                    out.write(f"{trg}\n")
                    
    elif strategy == 3:
        table_dic={}
        # Store the root distance of CDS and genome hits for each TRG.
        with open(table_path, 'r') as table:
            next(table)
            for line in table:
                line=line.strip().split('\t')
                
                focal=line[0] # name of the focal
                trg=line[1] # name of the TRG
                neighbor=line[2] # the smallest, the farest from the focal
                
                if trg not in table_dic:
                    table_dic[trg]=[]
                if neighbor != focal:
                    table_dic[trg].append(neighbor)
                    
        with open(out_path, 'w') as out:
            for trg in table_dic:
                # If there is not match outside of the focal (strict orphan)
                if table_dic[trg] == []:
                    trg=trg.rsplit('_elongated', 1)[0]
                    out.write(f"{trg}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", help="The path to the best hits table")
    parser.add_argument("--strategy", type=int, default=1, help="The strategy to apply")
    parser.add_argument("--synteny", default="False", help="Whether or not to check synteny for non-coding hits")
    parser.add_argument("--out", help="The output file")
    args = parser.parse_args()
    
    if args.synteny == True or args.synteny.upper() == "TRUE" :
        args.synteny = True
        
    filter_table(args.table, args.strategy, args.out, args.synteny)
