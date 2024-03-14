#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The idea behind this script is to identify alignments made of several HSP due to frameshifts in the subject.
When BLAST considers two or more HSP to be part of a same alignment, it recomputes a new evalue with sum-statistics.
all the different HSP of the same alignment thus share the same evalue but keep their own bitscore.

The tblastn command has been configured with -outfmt "7 std qlen qcovhsp sframe"
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, query length, % query coverage per hsp, sbjct frame
"""

import argparse
import sys
import math

# Define the command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help='Input tblastn format 7 file')
parser.add_argument('-o', '--outputfile', help='Output tblastn format 7 file (default: stdout)')

# Parse the command-line arguments
args = parser.parse_args()

# Get the column names
with open(args.inputfile, 'r') as inp:
        colnames = None
        for line in inp:
             if line.startswith('# Fields: '):
                 header = line.strip()[len('# Fields: '):]
                 colnames = header.split(', ')
                 break
        # print(f"Column names: {colnames}")

def sign(num):
    return int(math.copysign(1, int(num)))

# Open the input and the output files
with open(args.inputfile, 'r') as inp, (open(args.outputfile, 'w') if args.outputfile else sys.stdout) as out:

    if colnames is None: # If there is no hit at all in the file.
        print(f"No hit found in {args.inputfile}. The output file will be empty.")
        # Write a 0 octet file
        out.write('')
    else :

        previous_line = (". " * len(colnames)).split()
        previous_line = dict(zip(colnames, previous_line))

        for line in inp:
        # for i, line in enumerate(inp):
        #     print(f"line {i}: {line}")
            if line.startswith('#'):
                out.write(line)
            else:
                # print(line)
                line = line.strip().split('\t')
                line = dict(zip(colnames, line))
                try:
                    line['q. start'] = int(line['q. start'])
                except:
                    sys.exit(f"Error : non q.start for line: {line}")
                line['q. end'] = int(line['q. end'])
                line['s. start'] = int(line['s. start'])
                line['s. end'] = int(line['s. end'])
                line['query length'] = int(line['query length'])

                # If the query and subject are the same as for the previous HSP (line), and the evalue is non-null and the same as the previous HSP, and with the same strand but with a different bitscore,
                # then is it assumed that the current HSP belongs to the same alignment as the previous one.
                if line['query acc.ver'] == previous_line['query acc.ver'] and line['subject acc.ver'] == previous_line['subject acc.ver'] and line['evalue'] != "0" and line['evalue'] == previous_line['evalue'] and sign(line['sbjct frame']) == sign(previous_line['sbjct frame']) and line['bit score'] != previous_line['bit score']:
                    line['% identity'] = "NaN"
                    line['alignment length'] = "NaN"
                    line['mismatches'] = "NaN"
                    line['gap opens'] = "NaN"
                    line['bit score'] = "Nan"
                    line['q. start'] = min(line['q. start'], previous_line['q. start'])
                    line['q. end'] = max(line['q. end'], previous_line['q. end'])
                    line['% query coverage per hsp'] = round(100*(line['q. end'] - line['q. start'] + 1)/line['query length'])
                                        
                    strand = sign(line['sbjct frame'])
                    # If the strand is positive
                    if strand == 1:
                        line['sbjct frame'] = 42
                        line['s. start'] = min(line['s. start'], previous_line['s. start'])
                        line['s. end'] = max(line['s. end'], previous_line['s. end'])
                    # If the strand is negative
                    elif strand == -1:
                        line['sbjct frame'] = -42
                        line['s. start'] = max(line['s. start'], previous_line['s. start'])
                        line['s. end'] = min(line['s. end'], previous_line['s. end'])
                    else:
                        raise ValueError(f"Invalid strand: {strand}")
                # If the current HSP does not belong to the same alignment as the previous one
                else:
                    if previous_line['query acc.ver'] != '.':
                        outline = '\t'.join(str(previous_line[colname]) for colname in previous_line)
                        out.write(f"{outline}\n")
                        
                line['subject acc.ver'] = '_'.join( [ line['subject acc.ver'], str(line['s. start']), str(line['s. end']) ] )
                previous_line = line
    
    if colnames is not None:
        # Write the last line
        outline = '\t'.join(str(previous_line[colname]) for colname in previous_line)
        out.write(f"{outline}\n")