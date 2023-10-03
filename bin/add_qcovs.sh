#!/bin/bash
# Author : Paul Roginski
# Contact : plm.roginsk@gmail.com (if someones see an error or improve the script)
# Version : 230913_1010


# This script adds a qcovs column to a BLAST or DIAMONDv2 output.
# It was designed for this format : "6/7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
# Obviously one could pipe 2. | 3. and merge 3. and 4. but I'm not sure it would have a good perf/clarity ratio ;)
# '$2 == "*"' is for DIAMONDv2 outputs where it seems to mean that no hits were found the the query




# 1. Turn BLAST output to "BED"
awk '
	BEGIN{FS=OFS="\t"}

	# For non-comment line :
	# query-subject pair as sequence
	# qstart and qend as start and end
	!/#/ && $2!= "*" { print $1"_vs_"$2, $7-1, $8 }

' $1 > ${1}.tmp.bed




# 2. Use Bedtools to merge contiguous intervals
bedtools sort -i ${1}.tmp.bed > ${1}.tmp.bed.sorted
bedtools merge -i ${1}.tmp.bed.sorted > ${1}.tmp.bed.merged




# 3. Build a two cols .tsv file with query_vs_subject as col1 and the length part of the query covered by the target as col2
awk '
	BEGIN{FS=OFS="\t"}

	{ sum[$1]+= ($3-$2) }

	END{ for(key in sum){print key, sum[key]} }

' ${1}.tmp.bed.merged > ${1}.tmp_qlencov.tsv




# 4. Display the original alignment file with and additional column
awk '
	BEGIN{FS=OFS="\t"}

	# Store the first file into an array
	FNR==NR { qcovslen[$1] = $2 }

	# Second file (alingment)
	# Just print comment lines
	FNR!=NR && /^#/
	# For query with no hits
	FNR!=NR && !/^#/ && $2 == "*" { print $0,"*" }
	# For the other lines :
	FNR!=NR && !/^#/ && $2 != "*" {

		pair = $1"_vs_"$2

		# Interestingly, BLAST seems to provide rounded qcovs EXECEPT for 99.5 <= qcovs < 100  where it returns 99.
		# I suppose this is because it only wants to display 100% when EVERY single AA of the query is covered.
#		qcovs = int( (100*qcovslen[pair]/$13) +0.5)
		qcovs = 100*qcovslen[pair]/$13

		if(qcovs >= 99.5 && qcovs < 100){
			qcovs = 99
		} else {
			qcovs = int(qcovs + 0.5)
		}

		print $0, qcovs

	}

' ${1}.tmp_qlencov.tsv $1




rm ${1}.tmp.bed ${1}.tmp.bed.sorted ${1}.tmp.bed.merged ${1}.tmp_qlencov.tsv
