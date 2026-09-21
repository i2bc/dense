#!/bin/bash

# This script modifies inplace a nucleotide FASTA file by removing all sequences that do not end with a stop when translated.
# It also returns an aa version of the FASTA.

# WARNING : The script will keep the sequences that have a stop codon in the first frame, followed by one partial codon (unlikely).

fna="${1}"
faa=$(echo $fna | sed -E "s/(.*)\..*/\1/").faa

# Get a translated version of the input nucleotide FASTA
fastatranslate $fna -F 1 | sed "/>/ s/ \\[translate(1)\\]//" > $faa

# linearize the aa FASTA
awk '
	# credits to https://gist.github.com/lindenb/2c0d4e11fd8a96d4c345
	/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;}
	     {printf("%s",$0);}
	END  {printf("\n");}
' $faa > $faa.tmp

# Filter out sequences that do not end with a stop codon
awk '

    />/ { header=$0 }
    !/>/ && /\*$/ { print header ; print $0 }

' $faa.tmp > $faa

grep ">" $faa | sed "s/>//" > valid_sequences.txt

faSomeRecords $fna valid_sequences.txt ${fna}.tmp

# linearize the nucleotide FASTA
awk '
	# credits to https://gist.github.com/lindenb/2c0d4e11fd8a96d4c345
	/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;}
	     {printf("%s",$0);}
	END  {printf("\n");}
' ${fna}.tmp > $fna

rm valid_sequences.txt ${fna}.tmp ${faa}.tmp