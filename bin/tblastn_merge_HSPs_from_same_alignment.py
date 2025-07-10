#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
tblastn_merge_HSPs_from_same_alignment.py
-----------------------------------------

Merge HSPs from tblastn output that belong to the same alignment (frameshift-aware).

When BLAST considers two or more HSPs to be part of the same alignment, it recomputes a new e-value with sum-statistics.
All HSPs of the same alignment thus share the same e-value but keep their own bitscore.

This script identifies such HSPs (by matching query, subject, e-value, strand, and different bitscore) and merges them into a single output line,
extending the coordinates and updating the output subject identifier to include the merged region.

Required fields: query acc.ver, subject acc.ver, evalue, sbjct frame, bit score, q. start, q. end, s. start, s. end
All other columns are optional and preserved if present.

Author: Paul Roginski, 2025
License: MIT
"""

import argparse
import sys
import math
import warnings

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge tblastn HSPs from the same alignment (frameshift-aware)."
    )
    parser.add_argument('inputfile', help='Input tblastn format 7 file')
    parser.add_argument('-o', '--outputfile', help='Output tblastn format 7 file (default: stdout)')
    return parser.parse_args()

def get_colnames(inputfile):
    """Extract column names from tblastn format 7 header, return as list."""
    with open(inputfile, 'r') as inp:
        for line in inp:
            if line.startswith('# Fields: '):
                header = line.strip()[len('# Fields: '):]
                return [c.strip() for c in header.split(',')]
    warnings.warn("No '# Fields: ' header found in input file.")
    sys.exit(0)

def sign(num):
    """Return the sign of a number as int (+1 or -1)."""
    return int(math.copysign(1, int(num)))

def merge_tblastn_hsps(inputfile, outputfile=None):
    """
    Merge HSPs from tblastn output that belong to the same alignment (frameshift-aware).
    Requires the following fields in the input: query acc.ver, subject acc.ver, evalue, sbjct frame, bit score, q. start, q. end, s. start, s. end.
    Writes output to outputfile or stdout. The input can have any column order or extra columns.
    """
    required_fields = [
        'query acc.ver', 'subject acc.ver', 'evalue', 'sbjct frame', 'bit score',
        'q. start', 'q. end', 's. start', 's. end'
    ]
    colnames = get_colnames(inputfile)
    missing = [f for f in required_fields if f not in colnames]
    if missing:
        raise ValueError(f"Input file is missing required fields: {missing}")

    with open(inputfile, 'r') as inp, (open(outputfile, 'w') if outputfile else sys.stdout) as out:
        previous_line = (". " * len(colnames)).split()
        previous_line = dict(zip(colnames, previous_line))

        for line in inp:
            if line.startswith('#'):
                out.write(line)
            else:
                line = line.strip().split('\t')
                line = dict(zip(colnames, line))
                # Convert numeric fields (skip if not present)
                for field in ['q. start', 'q. end', 's. start', 's. end']:
                    if field in line:
                        line[field] = int(line[field])

                # Determine if this HSP can be merged with the previous one
                can_merge = (
                    previous_line['query acc.ver'] != '.' and
                    line['query acc.ver'] == previous_line['query acc.ver'] and
                    line['subject acc.ver'] == previous_line['subject acc.ver'] and
                    line['evalue'] != "0" and
                    line['evalue'] == previous_line['evalue'] and
                    sign(line['sbjct frame']) == sign(previous_line['sbjct frame']) and
                    line['bit score'] != previous_line['bit score']
                )

                if can_merge:
                    # Merge HSPs: update coordinates, set merged fields to NaN if present
                    for f in ['% identity', 'alignment length', 'mismatches', 'gap opens']:
                        if f in line: line[f] = "NaN"
                    line['bit score'] = "NaN"
                    line['q. start'] = min(line['q. start'], previous_line['q. start'])
                    line['q. end'] = max(line['q. end'], previous_line['q. end'])
                    # Only update % query coverage per hsp if present and query length is present
                    if '% query coverage per hsp' in line and 'query length' in line:
                        try:
                            qlen = int(line['query length'])
                            line['% query coverage per hsp'] = round(100 * (line['q. end'] - line['q. start'] + 1) / qlen)
                        except Exception:
                            line['% query coverage per hsp'] = "NaN"

                    strand = sign(line['sbjct frame'])
                    if strand == 1:
                        line['sbjct frame'] = 42
                        line['s. start'] = min(line['s. start'], previous_line['s. start'])
                        line['s. end'] = max(line['s. end'], previous_line['s. end'])
                    elif strand == -1:
                        line['sbjct frame'] = -42
                        line['s. start'] = max(line['s. start'], previous_line['s. start'])
                        line['s. end'] = min(line['s. end'], previous_line['s. end'])
                    else:
                        raise ValueError(f"Invalid strand: {strand}")
                    # Do not output previous line yet, keep merging
                else:
                    # Output previous line if it is not the dummy
                    if previous_line['query acc.ver'] != '.':
                        prev_out = previous_line.copy()
                        prev_out['subject acc.ver'] = '_'.join([
                            previous_line['subject acc.ver'],
                            str(previous_line['s. start']),
                            str(previous_line['s. end'])
                        ])
                        outline = '\t'.join(str(prev_out[colname]) for colname in colnames)
                        out.write(f"{outline}\n")

                previous_line = line

        # Write the last line
        if previous_line['query acc.ver'] != '.':
            prev_out = previous_line.copy()
            prev_out['subject acc.ver'] = '_'.join([
                previous_line['subject acc.ver'],
                str(previous_line['s. start']),
                str(previous_line['s. end'])
            ])
            outline = '\t'.join(str(prev_out[colname]) for colname in colnames)
            out.write(f"{outline}\n")

def main():
    args = parse_args()
    merge_tblastn_hsps(args.inputfile, args.outputfile)

if __name__ == "__main__":
    main()