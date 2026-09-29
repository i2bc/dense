#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check if the orthologous non-coding match for a de novo gene candidate has a disrupted integrity.

This is a conservative approach to verify the non-coding status of the dna sequence, and make sure
it isn't a gene that wouldn't have been annotated.

Inputs:
    - denovo_seq_aa = protein sequence of the de novo candidate
    - extended_match_seq_nt = dna sequence of the orthologous non-coding match in the outgroup, extended 
    99 nt (or the amount needed, multiple of 3) on each side.
    - integrity_threshold = threshold for the integrity check (default: 0.8). If the query coverage 
    of the alignment is below this threshold, the integrity is considered disrupted.


Author: Eliott Tempez, 2026
License: MIT
"""


import os
import tempfile
import argparse
import subprocess
import pandas as pd
from Bio.Seq import Seq
import psa


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check if the orthologous non-coding match for a de novo gene candidate has a disrupted integrity."
    )
    parser.add_argument('denovo_aa', help='Protein sequence of the de novo candidate')
    parser.add_argument('extended_nc_nt', help='DNA sequence of the extended orthologous non-coding match in the outgroup')
    parser.add_argument('--integrity_threshold', type=float, default=0.8, help='Threshold for the integrity check (default: 0.8)')
    return parser.parse_args()


def tblastn_from_sequences(query_sequence_aa, subject_sequence_nt):
    """Blast with tblastn a query sequence against a subject sequence, and return the blast result as a dataframe."""
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as query_file:
        query_file.write(f">query_aa\n{query_sequence_aa}\n")
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as subject_file:
        subject_file.write(f">subject_nt\n{subject_sequence_nt}\n")
    query_fasta_file = query_file.name
    subject_fasta_file = subject_file.name
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.tsv') as temp_output:
        output_file = temp_output.name
    
    out_command = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp sframe"
    command = [
            "tblastn",
            "-query", query_fasta_file,
            "-subject", subject_fasta_file,
            "-outfmt", out_command,
            "-out", output_file,
            "-evalue", "1e-3"
        ]
    subprocess.run(command, check=True, stdout=None, stderr=None)
    
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                "qlen", "qcovhsp", "sframe"]
    if os.path.getsize(output_file) == 0:
        return pd.DataFrame(columns=columns)
    blast_results = pd.read_csv(output_file, sep="\t", header=None, names=columns, comment="#")
    
    os.remove(query_fasta_file)
    os.remove(subject_fasta_file)
    os.remove(output_file)
    
    return blast_results


def local_alignment(qseq, sseq):
    """Perform a Smith-Waterman (water) alignment between two protein sequences and return the aligned sequences."""
    aln = psa.water(moltype="prot", qseq=qseq, sseq=sseq)
    return aln


def get_qcov_from_alignment(aln, qseq):
    """Extract the query coverage from the local alignment obtained with psa."""
    aln_length = aln.length
    query_length = len(qseq)
    qcov = aln_length / query_length
    return qcov


def integrity_search(denovo_seq_aa, extended_match_seq_nt, integrity_threshold):
    ## 1 - Frameshift
    blast_result = tblastn_from_sequences(denovo_seq_aa, extended_match_seq_nt)
    blast_result = blast_result[blast_result["sframe"] > 0].reset_index(drop=True)
    if blast_result.empty:
        return "no_match"
    if len(blast_result.index) > 1:
        return "frameshift"
    
    frame_with_best_qcov = blast_result.loc[blast_result["qcovhsp"].idxmax()]["sframe"]
        
    ## 2 - Truncated
    extended_match_seq_aa = Seq(extended_match_seq_nt[frame_with_best_qcov - 1:]).translate(table=11)
    aln = local_alignment(denovo_seq_aa, extended_match_seq_aa)
    qcov = get_qcov_from_alignment(aln, denovo_seq_aa)
    if qcov < integrity_threshold:
        return "truncated"
    
    ## 3 - Internal stop
    aln_subject_seq = aln.sseq
    denovo_length = len(denovo_seq_aa)
    fragmented_subject_seq = aln_subject_seq.split("*")
    biggest_chunk_length = max([len(chunk) for chunk in fragmented_subject_seq])
    if biggest_chunk_length / denovo_length < integrity_threshold:
        return "internal_stop"
    
    ## 4 - No disruption
    return "intact"


def main():
    args = parse_args()
    integrity_status = integrity_search(args.denovo_aa, args.extended_nc_nt, args.integrity_threshold)
    print(integrity_status)


if __name__ == "__main__":
    main()
