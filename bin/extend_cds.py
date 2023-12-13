#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pybedtools.bedtool import BedTool,Interval
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess


def extend_cds(out, gfasta, fai, gff, cds, size=99, mode='simple'):
    """
    If you want the flanking nucleotides of a genomic feature you should probably use bedtools flank.
    However, since CDS can be represented by multiple genomique features in a GFF file, you can use this program to first find its real borders before using bedtools flank.
    # This script get the flanking n nucleotides at each side of a given CDS (the n upstream and n downstream nucleotides with respect to strand).
    # Optionnaly, it then multitranslate it so that the final FASTA contains for each CDS :
    	# n_nucl_translated_in_frame_0____translated_CDS____n_nucl_translated_in_frame_0
    	# n_nucl_translated_in_frame_0____translated_CDS____n_nucl_translated_in_frame_1
    	# n_nucl_translated_in_frame_0____translated_CDS____n_nucl_translated_in_frame_2
    	# n_nucl_translated_in_frame_1____translated_CDS____n_nucl_translated_in_frame_0
    	# n_nucl_translated_in_frame_1____translated_CDS____n_nucl_translated_in_frame_1
    	# n_nucl_translated_in_frame_1____translated_CDS____n_nucl_translated_in_frame_2
    	# n_nucl_translated_in_frame_2____translated_CDS____n_nucl_translated_in_frame_0
    	# n_nucl_translated_in_frame_2____translated_CDS____n_nucl_translated_in_frame_1
    	# n_nucl_translated_in_frame_2____translated_CDS____n_nucl_translated_in_frame_2
    
    # It takes as input (all mandatory) :
    # - a genomic fasta file (gfasta),
    # - the corresponding GFF3 file (gff),
    # - the genome index file (created if necessary) containning the size of each sequence ("chromosome") from the genomic fasta (fai),
    # - a nucleotide fasta to elongated (header must match an 'mRNA' feature's 'ID' (col9)) (cds),
    # - the number of nucleotides to retrieve on each size of the CDS (size),
    # - whether or not to perform a multiframe translation (mode)
    
    # WARNING #
    # The GFF3 must be valid with start (col3) < stop (col4).
    # The element to elongate must be the 'Parent' (col9) of 'CDS' (col3) features.
    # CDS features must have a strand (col7).
    # Multiple CDS features belonging to the same 'Parent' must have the same strand.
    # If this is not the case (rare), the first (closest to 5') CDS feature will serve to define the unique strand, and thus the output will be nonsense. 

    # To get a fai
    # samtools faidx GRCh38.fa
    # bedtools flank -i my.bed -g GRCh38.fa.fai
    """

    if size%3 != 0:
        print(f"WARNING : {size} is not a multiple of three. If you translate the prefix you will not be in the main sequence's frame.")
    
    # From the cds file, build a dictionnary with CDS' names as keys.
    parent_dic={}
    with open(cds, "r") as query:
        for record in SeqIO.parse(query, "fasta"):
            parent_dic[record.name]={'seq':record.seq}
#    parents=list(parent_dic.keys())

    #Quickly reduce the gff file (nothing goes faster than grep for that)
    gff_CDS=subprocess.check_output(["grep", "CDS", gff])
    gff_CDS = gff_CDS.decode("ascii",errors="ignore")

    # In the gff, look for features of type "CDS" whose "Parent" attribute is in the list of CDS to elongate.
    pattern = r"Parent=([^;]+)"
#    with open(gff, "r") as gff_file:
#        for line in gff_file:
    for line in gff_CDS.splitlines():
        if not line.startswith("#"):  # Skip header lines
            entry = line.strip().split("\t")
            if entry[2] == "CDS":
                match = re.search(pattern, entry[8])
                if match and match.group(1) in parent_dic:
                    parent=match.group(1)
                    if not 'strand' in parent_dic[parent]: #If the parent strand is not set, use the feature's strand to define it. This implies that in the unlikely case where a CDS is made of features of different strand, the first strand encountered will be use.
                        parent_dic[parent]['strand']=entry[6]
                    if not 'chrom' in parent_dic[parent]: #If the parent strand is not set, use the feature's strand to define it. This implies that in the unlikely case where a CDS is made of features of different strand, the first strand encountered will be use.
                        parent_dic[parent]['chrom']=entry[0]
                    if not 'start' in parent_dic[parent] or parent_dic[parent]['start'] > int(entry[3]): # Look for the minimum start position among the 'CDS' features sharing the parent
                        parent_dic[parent]['start']=int(entry[3])
                    if not 'end' in parent_dic[parent] or parent_dic[parent]['end'] < int(entry[4]): # Look for the minimum start position among the 'CDS' features sharing the parent
                        parent_dic[parent]['end']=int(entry[4])
    
    #Build a bedtool object from the collected features to use bedtools flank and thus retrieve the n flanking nucleotides.
    bed_entries = []
    for parent in parent_dic:
        entry = parent_dic[parent]
        chrom = entry['chrom']
        start = entry['start']-1
        end = entry['end']
        strand = entry['strand']
        bed_entry = Interval(chrom=chrom, start=start, end=end, name=parent, score=".", strand=strand)
        bed_entries.append(bed_entry)
    
    bed = BedTool(bed_entries)
    flanked = bed.flank(g=fai, b=size)
    flanks = flanked.sequence(fi = gfasta, s=True, nameOnly=True) #Make a tmp fasta out of it (equivalent to bedtools getfasta)
    
    # Define which flank is prefix or suffix according to the strand
    with open(flanks.seqfn, "r") as flanks_fna:
        for record in SeqIO.parse(flanks_fna, "fasta"):
            parent = record.name[:-3]
            sequence = record.seq
            if not 'prefix' in parent_dic[parent] and not 'suffix' in parent_dic[parent]:
                if parent_dic[parent]['strand'] == "+":
                    parent_dic[parent]['prefix'] = sequence
                elif parent_dic[parent]['strand'] == "-":
                    parent_dic[parent]['suffix'] = sequence
            else :
                if 'prefix' in parent_dic[parent]:
                    parent_dic[parent]['suffix'] = sequence
                elif 'suffix' in parent_dic[parent]:
                    parent_dic[parent]['prefix'] = sequence
                
      
    def translate_with_frame(seq, frame):
        
        if frame not in range(3):
            raise ValueError("Frame must be between 0 and 2")
    
        # Pad the sequence with 'N's to ensure it's a multiple of three
        padded_sequence = seq[frame:]+'N'*frame
        
        # Translate the padded sequence in frame 0 (no shift)
        translated_sequence = padded_sequence.translate()
    
        # Remove the additional 'X' from the translated sequence
        if frame in [1,2]:
            translated_sequence = translated_sequence[:-1]
    
        return translated_sequence
    
    #Write the output file according to the required mode.
    with open(out, 'w') as output: 
        if mode == "nucl" :
            for parent in parent_dic :
                #Return a fasta with prefix->CDS->suffix translated in frame 0
                output.write(f">{parent}_elongated {size} nucleotides\n")
                if "prefix" in parent_dic[parent]:
                    output.write(f"{parent_dic[parent]['prefix']}")
                output.write(f"{parent_dic[parent]['seq']}")
                if "suffix" in parent_dic[parent]:
                    output.write(f"{parent_dic[parent]['suffix']}")
                output.write(f"\n")
                # output.write(f"{parent_dic[parent]['prefix']}{parent_dic[parent]['seq']}{parent_dic[parent]['suffix']}\n")
        
        elif mode == "simple" :
            for parent in parent_dic :
                #Return a fasta with prefix->CDS->suffix translated in frame 0
                output.write(f">{parent}_elongated {size} nucleotides\n")
                # output.write(f"{parent_dic[parent]['prefix'].translate()}{parent_dic[parent]['seq'].translate()}{parent_dic[parent]['suffix'].translate()}\n")
                if "prefix" in parent_dic[parent]:
                    output.write(f"{parent_dic[parent]['prefix'].translate()}")
                output.write(f"{parent_dic[parent]['seq'].translate()}")
                if "suffix" in parent_dic[parent]:
                    output.write(f"{parent_dic[parent]['suffix'].translate()}")
                output.write(f"\n")

        elif mode == "multi" :
            for parent in parent_dic :
                #Return a fasta with 9 entries per original CDS (3 frames possible for prefix x 3 frames possible for suffix)
                for pf in range(3):
                    for sf in range(3):
                        output.write(f">{parent}_elongated {size} nucleotides pf:{pf} sp:{sf}\n")
                        if "prefix" in parent_dic[parent]:
                            prefix_trans=translate_with_frame(parent_dic[parent]['prefix'], pf)
                            output.write(f"{prefix_trans}")
                        output.write(f"{parent_dic[parent]['seq'].translate()}")
                        if "suffix" in parent_dic[parent]:
                            suffix_trans=translate_with_frame(parent_dic[parent]['suffix'], sf)
                            output.write(f"{suffix_trans}")
                        output.write(f"\n")
                        # output.write(f"{prefix_trans}{parent_dic[parent]['seq'].translate()}{suffix_trans}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gfasta", help="The path to the genomic fasta file")
    parser.add_argument("--fai", help="The path to the genomic fasta index file")
    parser.add_argument("--gff", help="The path to the corresponding GFF3 file")
    parser.add_argument("--cds", help="The path to the fasta file containing CDS sequences")
    parser.add_argument("--size", type=int, help="The number of nucleotides to add on each side of the CDS")
    parser.add_argument("--mode", help="nucl : Print the nucleotide elongated CDS\nsimple : Print the amino acid translation of the elongated CDS\nmulti : Return 9 entries per CDS (3x3 frames for prefix and suffix)")
    parser.add_argument("--out", help="output fasta")
    args = parser.parse_args()
    
#    HOME = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/" 
#    query_path = HOME+"Eukaryotes/MODELS/SCER/Scer_NCBI_TRG_10.fna"
#    fasta_path = HOME+"Eukaryotes/MODELS/SCER/GENOMES/Scer_NCBI.fna"
#    fai_path = HOME+"Eukaryotes/MODELS/SCER/GENOMES/Scer_NCBI.fna.fai"
#    gff_filename = HOME+"Eukaryotes/MODELS/SCER/GENOMES/Scer_NCBI.gff"
#    
#    HOME = "/home/paul.roginski/Bureau/"
#    query_path = HOME+"Scer_NCBI_TRG_10.fna"
#    fasta_path = HOME+"Scer_NCBI.fna"
#    fai_path = HOME+"Scer_NCBI.fna.fai"
#    gff_filename = HOME+"Scer_NCBI.gff"
#    mode="multiple"
#    size=99

    extend_cds(args.out, args.gfasta, args.fai, args.gff, args.cds, args.size, args.mode)