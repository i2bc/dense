taxdump=$1
taxid=$2
TRG_rank=$3

rank=""
counter=0

while [[ $rank != $TRG_rank ]] && [[ $rank != "no rank" ]] && [[ $counter -lt 100 ]]
do
	old_taxid=$taxid
	rank=$( grep -P "^${taxid}\t" ${taxdump}/nodes.dmp | sed -E "s/([0-9]+)\t\|\t([0-9]+)\t\|\t([^\t]+)\t.*/\3/")
	taxid=$(grep -P "^${taxid}\t" ${taxdump}/nodes.dmp | sed -E "s/([0-9]+)\t\|\t([0-9]+)\t\|\t([^\t]+)\t.*/\2/")
	counter=$((( $counter + 1 )))
#	echo $old_taxid $rank $counter
done

if [[ $rank != "no rank" ]] && [[ ! $counter -eq 100 ]]
then
	grep -P "^${old_taxid}\t" ${taxdump}/fullnamelineage.dmp | sed -E "s/([0-9]+)\t\|\t([^\t]+)\t\|.*/\2/"
else
	echo "${TRG_rank} could not be found in ${taxid} phylogeny."
	exit 1
fi
