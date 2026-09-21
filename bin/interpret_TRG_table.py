#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:55:42 2023

@author: paul.roginski
"""

import argparse
import csv

def table_summary(table_path, out_path):
    """
    This function builds a summary of the main TRG_table.
    """

    table_dic={} # to store the name and the distance of each TRG's hits
    with open(table_path, 'r') as table:
        next(table)
                
        for line in table:
            line=line.strip().split('\t')
            
            trg=line[1] # name of the TRG
            name=line[2] # name of the 'neighbor' species
            hit_type=line[4] # CDS or genome
            root_distance=float(line[6]) # the smallest, the farest from the focal
            
            if trg not in table_dic:
                table_dic[trg]={
                        'CDS_hits_root_distance' : [],
                        'CDS_hits_name' : [],
                        'genome_hits_root_distance' : [],
                        'genome_hits_name' : [],
                        'genome_hits_with_synteny_root_distance' : []
                        }
                
            if hit_type == "CDS":
                table_dic[trg]['CDS_hits_root_distance'].append(root_distance)
                table_dic[trg]['CDS_hits_name'].append(name)
                
            elif hit_type == "genome":
                table_dic[trg]['genome_hits_root_distance'].append(root_distance)
                table_dic[trg]['genome_hits_name'].append(name)

                if line[7] == "True" :# Boolean
                    table_dic[trg]['genome_hits_with_synteny_root_distance'].append(root_distance)
                    
                    
    # Get the focal species' name
    with open(table_path, 'r') as table:
        next(table)
        focal = table.readline().strip().split('\t')[0]  
    
    
    # Write a summary TSV file
    with open(out_path, 'w') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(["TRG_CDS", "isOrphan", "hasNcHit", "hasNcHitSynt", "hasNcHitOutgroup", "hasNcHitOutgroupSynt"]) # header
        
        for trg in table_dic :
            sub=table_dic[trg]
            
            if sub['CDS_hits_name'] == [focal] :
                isOrphan = True
            else :
                isOrphan = False
            
            nc_hit_names = list(set(sub['genome_hits_name']) - set(sub['CDS_hits_name']))
            if nc_hit_names != []:
                hasNcHit = True
            else :
                hasNcHit = False
                
            if sub['genome_hits_with_synteny_root_distance'] != []:
                hasNcHitSynt = True
            else :
                hasNcHitSynt = False
            
            farest_CDS=float('inf')
            farest_genome=float('inf')
            if len(sub['CDS_hits_root_distance']) > 0 :
                farest_CDS=min(sub['CDS_hits_root_distance'])
            if len(sub['genome_hits_root_distance']) > 0 :
                farest_genome=min(sub['genome_hits_root_distance'])
                    
            # If the farest genome hits is farer than the farest CDS hit
            if farest_genome < farest_CDS:
                hasNcHitOutgroup = True
            else :
                hasNcHitOutgroup = False
                
            farest_genome_synt=float('inf')
            if len(sub['genome_hits_with_synteny_root_distance']) > 0 :
                farest_genome_synt=min(sub['genome_hits_with_synteny_root_distance'])
                
            # If the farest genome hits in synteny is farer than the farest CDS hit
            if farest_genome_synt < farest_CDS:
                hasNcHitOutgroupSynt = True
            else :
                hasNcHitOutgroupSynt = False

            trg_cds = trg.rsplit('_elongated', 1)[0]

            writer.writerow([trg_cds, isOrphan, hasNcHit, hasNcHitSynt, hasNcHitOutgroup, hasNcHitOutgroupSynt])
            



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", help="The path to the best hits table")
    parser.add_argument("--out", help="The output file")
    args = parser.parse_args()
    
    table_summary(args.table, args.out)
