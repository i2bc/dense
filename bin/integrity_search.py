#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check if the orthologous non-coding match for a de novo gene candidate has a disrupted integrity.

This is a conservative approach to verify the non-coding status of the dna sequence, and make sure
it isn't a gene that wouldn't have been annotated.

Inputs:
    - denovo_aa_fasta      : FASTA file with the protein sequence of the de novo candidate
    - extended_nc_nt_fasta : FASTA file with the DNA sequence of the extended orthologous
                             non-coding match in the outgroup (extended 99 nt on each side).
    - --integrity_threshold : threshold for the integrity check (default: 0.8). If the query
      coverage of the alignment is below this threshold, the integrity is considered disrupted.

Author: Eliott Tempez, 2026
License: MIT
"""

import os
import tempfile
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check if the orthologous non-coding match for a de novo gene candidate has a disrupted integrity."
    )
    parser.add_argument('denovo_aa_fasta',
                        help='FASTA file with the protein sequence of the de novo candidate')
    parser.add_argument('extended_nc_nt_fasta',
                        help='FASTA file with the DNA sequence of the extended orthologous non-coding match')
    parser.add_argument('--integrity_threshold', type=float, default=0.8,
                        help='Threshold for the integrity check (default: 0.8)')
    return parser.parse_args()


def read_first_sequence(fasta_file):
    """Read and return the first sequence from a FASTA file as a string."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        return str(record.seq)
    raise ValueError(f"No sequence found in {fasta_file}")


def tblastn_from_files(query_fasta, subject_fasta):
    """Run tblastn with query (AA) and subject (NT) FASTA files; return results as a DataFrame."""
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.tsv') as tmp:
        output_file = tmp.name

    out_fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp sframe"
    command = [
        "tblastn",
        "-query",   query_fasta,
        "-subject", subject_fasta,
        "-outfmt",  out_fmt,
        "-out",     output_file,
        "-evalue",  "1e-3"
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore",
               "qlen", "qcovhsp", "sframe"]
    if os.path.getsize(output_file) == 0:
        os.remove(output_file)
        return pd.DataFrame(columns=columns)
    df = pd.read_csv(output_file, sep="\t", header=None, names=columns, comment="#")
    os.remove(output_file)
    return df


def local_alignment_qcov(qseq, sseq):
    """
    Perform a Smith-Waterman (local) alignment between two protein sequences
    using Bio.Align.PairwiseAligner, and return:
      - qcov  : query coverage (aligned_length / query_length)
      - saln  : aligned subject sequence (with gaps removed, stops kept)
    """
    from Bio.Align import substitution_matrices as sm

    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = sm.load("BLOSUM62")
    aligner.open_gap_score   = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(qseq, sseq)
    try:
        aln = next(iter(alignments))
    except StopIteration:
        return 0.0, ""

    aln_length = aln.shape[1]
    qcov = aln_length / len(qseq)

    aln_str = str(aln)
    lines = aln_str.strip().split('\n')
    saln = lines[0].split()[2] if len(lines) >= 3 else ""
    return qcov, saln


def integrity_search(denovo_aa_fasta, extended_nc_nt_fasta, integrity_threshold):
    denovo_seq_aa        = read_first_sequence(denovo_aa_fasta)
    extended_match_seq_nt = read_first_sequence(extended_nc_nt_fasta)

    ## 1 - Frameshift
    blast_result = tblastn_from_files(denovo_aa_fasta, extended_nc_nt_fasta)
    blast_result = blast_result[blast_result["sframe"] > 0].reset_index(drop=True)
    if blast_result.empty:
        return "no_match"
    if len(blast_result.index) > 1:
        return "frameshift"

    frame_with_best_qcov = int(blast_result.loc[blast_result["qcovhsp"].idxmax()]["sframe"])

    ## 2 - Truncated
    extended_match_seq_aa = str(Seq(extended_match_seq_nt[frame_with_best_qcov - 1:]).translate(table=11))
    qcov, saln = local_alignment_qcov(denovo_seq_aa, extended_match_seq_aa)
    if qcov < integrity_threshold:
        return "truncated"

    ## 3 - Internal stop
    denovo_length = len(denovo_seq_aa)
    fragmented = saln.replace('-', '').split("*")
    biggest_chunk_length = max((len(chunk) for chunk in fragmented), default=0)
    if denovo_length == 0 or biggest_chunk_length / denovo_length < integrity_threshold:
        return "internal_stop"

    ## 4 - No disruption
    return "intact"


def main():
    args = parse_args()
    result = integrity_search(args.denovo_aa_fasta, args.extended_nc_nt_fasta,
                              args.integrity_threshold)
    print(result)


if __name__ == "__main__":
    main()
