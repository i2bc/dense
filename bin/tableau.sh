#!/bin/bash

# The tblastn command has been configured with -outfmt "7 std qlen qcov(hsp/s) sframe".

# Get bash named parameters (adapted from https://www.brianchildress.co/named-parameters-in-bash/).
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
  shift
done


#touch ${in}.best_hits.tsv
touch TRG_blast_${subject_genome}_best_hits.tsv


awk -v query_genome="${query_genome}" -v subject_genome="${subject_genome}" -v type="${type}" '


	BEGIN{ FS=OFS="\t" }


	# For non-comment lines with a query coverage by subject  >= 50% (the evalue is assumed to have been filtered within the blast command),
	$1 !~ /#/ && $14 >= 50 {

		# Get the query name but without the specified frames.
		$1=gensub(/(.*_elongated)_F[0-9]_[0-9]/,"\\1","g",$1)

		# If this is the first time the programme meet a significant aligment for this query and this subject, add them to the "hits.tsv" file.
		if ( !($1$2 in data)){
			# Save the pair.
			data[$1$2]=1
			#print query_genome,$1,subject_genome,$2,type >> FILENAME".hits.tsv"
			print query_genome,$1,subject_genome,$2,type >> "TRG_blast_"subject_genome"_hits.tsv"
		}

		# If this alignment has the best evalue seen for the current query,
			if ( !($1 in evalue) || $11 < evalue[$1] ){
			# Record the evalue,
			evalue[$1]=$11
			# Record the subject.
			best[$1]=$2
		}
	}


	# At the end, for each query with a best-subject, add the query-subject pair to the "best_hits.tsv" file.
	END {
		for (query in best){
			#print query_genome,query,subject_genome,best[query],type >> FILENAME".best_hits.tsv"
			print query_genome,query,subject_genome,best[query],type >> "TRG_blast_"subject_genome"_best_hits.tsv"
		}
	}


' $in