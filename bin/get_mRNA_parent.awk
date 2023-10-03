BEGIN{FS=OFS="\t"}

# For lines that are :
# non-comment,
# of type "mRNA",
# with both "ID" and "Parent" in their attributes,
$0 !~ /^#/ && $3 ~ /^mRNA$/ && $9 ~ /ID=/ && $9 ~ /Parent=/ {

	ID=gensub(/.*ID=([^;]+).*/,"\\1","g",$9)

	# If this is the first time the mRNA is met, print the mRNA-gene pair and save the ID.
	if (!( ID in data )) {

	        mRNA_parent=gensub(/.*Parent=([^;]+).*/,"\\1","g",$9)

		print ID, mRNA_parent

		data[ID] = 1
	}

}
