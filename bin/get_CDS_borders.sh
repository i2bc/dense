#!/bin/bash

# This script generates two sequences for each CDS : the n upstream and n downstream nucleotides (with respect to strand).
# It then multitranslate it so that the final FASTA contains for each CDS :
	# 100_nucl_translated_in_frame_0____translated_CDS____100_nucl_translated_in_frame_0
	# 100_nucl_translated_in_frame_0____translated_CDS____100_nucl_translated_in_frame_1
	# 100_nucl_translated_in_frame_0____translated_CDS____100_nucl_translated_in_frame_2
	# 100_nucl_translated_in_frame_1____translated_CDS____100_nucl_translated_in_frame_0
	# 100_nucl_translated_in_frame_1____translated_CDS____100_nucl_translated_in_frame_1
	# 100_nucl_translated_in_frame_1____translated_CDS____100_nucl_translated_in_frame_2
	# 100_nucl_translated_in_frame_2____translated_CDS____100_nucl_translated_in_frame_0
	# 100_nucl_translated_in_frame_2____translated_CDS____100_nucl_translated_in_frame_1
	# 100_nucl_translated_in_frame_2____translated_CDS____100_nucl_translated_in_frame_2

# It takes as input (all mandatory) :
# - a genomic fasta file (gfasta)
# - the corresponding GFF file (gff)
# - the name of the genome file (created if necessary) containning the size of each sequences from the genomic fasta (genome)
# - the number of nucleotides to retrieve on each size of the CDS (size)
# - whether of not to perform a multiframe translation (multiframe)


# Get bash named parameters (adapted from https://www.brianchildress.co/named-parameters-in-bash/)
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
  shift
done


# Radical name that will be used for output files
name=$(basename $gfasta .fna)
pref_suf="${name}_prefix_n_suffix"


# Given a size,
# Among features of the gff having the same mRNA as parent,
# determine the minimum and the maximum coordinates (start and end of the CDS),
# store the seq, the source, and the strand to generate a factice GFF file.
# With respect to the strand, define a prefix (nucleotides before the CDS) and a suffix (after) of length "size"
# (limited by the bounds of the sequence (chromosome)).
awk -v size="${size}" '

function max(a, b){
   if (a > b)
   return a
   return b
}
function min(a, b){
   if (a < b)
   return a
   return b
}

BEGIN{
FS=OFS="\t"
# Minimum size of the prefix and suffix (e.g. 3 to have at least of translatable triplet).
min_size=3
size=max(size,min_size)
}

{
# For each line if the GFF file with "CDS" as feature type and with a "Parent=" attribute in the 9th field :
if( $3 == "CDS" && $9 ~ "Parent=" ) {

	parent=gensub(/.*Parent=([^;]+);.*/,"\\1","g",$9)

	if(parent in parents){
		if ($4 < start[parent]){start[parent]=$4}
		if ($5 > end[parent]){end[parent]=$5}
	}
	else {
		parents[parent]=1
		start[parent]=$4
		end[parent]=$5
		seq[parent]=$1
		source[parent]=$2
		strand[parent]=$7
	}

}
# For each line of the second file (= genome file):
if(FNR!=NR){ chr_length[$1] = $2 }
}

END {
for(parent in parents){

	# start position of the upstream segment to define
	up_start = start[parent]-size
        # end position of the upstream segment to define
	up_end = start[parent]-1
        # start position of the downstream segment to define
	down_start = end[parent]+1
        # end position of the downstream segment to define
	down_end = end[parent]+size

	# Print a "fake" gff line for the prefix and the suffix
	if(strand[parent] == "+"){
		if(up_end >= min_size){ 				print seq[parent],source[parent],parent"_prefix", max(1, up_start) , up_end ,".",strand[parent],".",parent"="start[parent]":"end[parent] }
		else{print "#"seq[parent],source[parent],parent"_prefix", max(1, up_start) , up_end ,".",strand[parent],".",parent"="start[parent]":"end[parent] "up_end="up_end ";min_size="min_size }
		if(down_start <= chr_length[seq[parent]] - min_size+1){ 	print seq[parent],source[parent],parent"_suffix", down_start , min(chr_length[seq[parent]], down_end) ,".",strand[parent],".",parent"="start[parent]":"end[parent] }
		else{print "#"seq[parent],source[parent],parent"_suffix", down_start , min(chr_length[seq[parent]], down_end) ,".",strand[parent],".",parent"="start[parent]":"end[parent] "down_start="down_start ";chr_length[seq[parent]]="chr_length[seq[parent]]}
	}
	else{
		if(down_start <= chr_length[seq[parent]] - min_size+1){	print seq[parent],source[parent],parent"_prefix", down_start , min(chr_length[seq[parent]], down_end) ,".",strand[parent],".",parent"="start[parent]":"end[parent]}
		else{print "#"seq[parent],source[parent],parent"_prefix", down_start , min(chr_length[seq[parent]], down_end) ,".",strand[parent],".",parent"="start[parent]":"end[parent] "down_start="down_start ";chr_length[seq[parent]]="chr_length[seq[parent]]}
	        if(up_end >= min_size){ 				print seq[parent],source[parent],parent"_suffix", max(1, up_start) , up_end ,".",strand[parent],".",parent"="start[parent]":"end[parent] }
	        else{print "#"seq[parent],source[parent],parent"_suffix", max(1, up_start) , up_end ,".",strand[parent],".",parent"="start[parent]":"end[parent] "up_end="up_end ";min_size="min_size}
	}

}
}' $gff $genome > ${pref_suf}.gff


# From the "fake" GFF file ${pref_suf}.gff, extract the prefix and the suffix sequences.
bedtools getfasta -s -nameOnly -fi $gfasta -bed ${pref_suf}.gff > ${pref_suf}.fna

# Translate the nucleotides FASTA with prefixes and suffixes into the three frames.
rm -f ${pref_suf}.faa
if [ $multiframe == "True" ]; then

	for frame in 0 1 2; do

		faTrans ${pref_suf}.fna ${pref_suf}_F${frame}.faa -offset=${frame}
		# Turn "Z" into "*", addd the frame information and linearize the FASTA
		sed "/^[^>]/ s/Z/*/g" ${pref_suf}_F${frame}.faa | \
		sed "s/>.*/&_F${frame}/" | \
		awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' >> ${pref_suf}.faa

		rm ${pref_suf}_F${frame}.faa

	done

else

	faTrans ${pref_suf}.fna ${pref_suf}.faa.tmp
	# Turn "Z" into "*", addd the frame information and linearize the FASTA
	sed "/^[^>]/ s/Z/*/g" ${pref_suf}.faa.tmp | \
	awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' > ${pref_suf}.faa

	rm ${pref_suf}.faa.tmp

fi



# Then according the the "--multiframe" option, build an elongated translated CDS FASTA.
echo "multiframe $multiframe"
if [ $multiframe == "True" ]; then
	awk '

	BEGIN{FS=OFS="\t"}

	{
	# For header lines of the first file :
	if(FNR==NR && $0 ~ /^>/){

		# Retrieve the CDS parent in e.g. >Sbay_4.469_prefix(-)_F0
		parent=gensub(/>(.*)_([a-z]+)\(.\)_F([0-9])/,"\\1","g",$0)
		side=gensub(/>(.*)_([a-z]+)\(.\)_F([0-9])/,"\\2","g",$0)
		frame=gensub(/>(.*)_([a-z]+)\(.\)_F([0-9])/,"\\3","g",$0)

		# Store the following line as the prefix/suffix for this parent
		getline
		if(side == "prefix"){ prefix[parent"_F"frame]=$0 }
		if(side == "suffix"){ suffix[parent"_F"frame]=$0 }


	}

	# For header lines of the second file :
	if(FNR!=NR ){

		# If the line is a header save it as header
		if($0 ~ /^>/){ header=$0 }
		# If the line is a seq :
		else {

			# Retrieve the CDS parent
			parent=gensub(/>/,"","g",header)

			# For the 9 possible combinaisons of prefix and suffix in the three frames :
			for (frame_pref=0; frame_pref<=2; frame_pref++){
				for (frame_suff=0; frame_suff<=2; frame_suff++){

					print header "_elongated_F"frame_pref"_"frame_suff
					print prefix[parent"_F"frame_pref] $0 suffix[parent"_F"frame_suff]

				}
			}
		}
	}
	}
	' ${pref_suf}.faa TRG.faa >> TRG_multielongated.faa
else
	awk '
	BEGIN{FS=OFS="\t"}
	{
	# For header lines of the first file :
	if(FNR==NR && $0 ~ /^>/){

		# Retrieve the CDS parent in e.g. >Sbay_4.469_prefix(-)
		parent=gensub(/>(.*)_([a-z]+)\(.\)/,"\\1","g",$0)
		side=gensub(/>(.*)_([a-z]+)\(.\)/,"\\2","g",$0)

		# Store the following line as the prefix/suffix for this parent
		getline
		if(side == "prefix"){ prefix[parent]=$0 }
		if(side == "suffix"){ suffix[parent]=$0 }


	}
	# For header lines of the second file :
	if(FNR!=NR ){

		# If the line is a header save it as header
		if($0 ~ /^>/){ header=$0 }
		# If the line is a seq :
		else {

			# Retrieve the CDS parent
			parent=gensub(/>/,"","g",header)

			print header "_elongated"
			print prefix[parent] $0 suffix[parent]

		}
	}
	}
	' ${pref_suf}.faa CDS.faa >> CDS_elongated.faa
fi



# Finally, also build an elongated CDS FASTA.
# It will be further used to make nucleotide aligments.
if [ $multiframe != "True" ]
then
	awk '
	BEGIN{FS=OFS="\t"}
	{
	## For header lines of the first file :
	if(FNR==NR && $0 ~ /^>/){

		# Retrieve the CDS parent in e.g. >Sbay_4.469_prefix(-)
		parent=gensub(/>(.*)_([a-z]+)\(.\)/,"\\1","g",$0)
		side=gensub(/>(.*)_([a-z]+)\(.\)/,"\\2","g",$0)

		# Store the following line as the prefix/suffix for this parent
		getline
		if(side == "prefix"){ prefix[parent]=$0 }
		if(side == "suffix"){ suffix[parent]=$0 }


	}
	# For header lines of the second file :
	if(FNR!=NR ){

		# If the line is a header save it as header
		if($0 ~ /^>/){ header=$0 }
		# If the line is a seq :
		else {

			# Retrieve the CDS parent
			parent=gensub(/>/,"","g",header)

			print header "_elongated"
			print prefix[parent] $0 suffix[parent]

		}
	}
	}
	' ${pref_suf}.fna CDS.fna >> CDS_elongated.fna
fi
