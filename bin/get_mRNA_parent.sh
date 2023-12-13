#!/bin/bash

# Look for an mRNA feature in the GFF
hasmRNA=$( awk -F"\t" 'BEGIN {out = "false"} $3 ~ /^mRNA$/ {out = "true"; exit 0} END {print out}' $1)

if [ "${hasmRNA}" == "true" ]
then

	# $1 has mRNA feature
	awk '
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
	' $1

else
    # $1 has no mRNA feature
	awk '
                BEGIN{FS=OFS="\t"}

                # For lines that are :
                # non-comment,
                # of type "gene",
                $0 !~ /^#/ && $3 ~ /^gene$/ {

                        ID=gensub(/.*ID=([^;]+).*/,"\\1","g",$9)

                        # If this is the first time the gene is met, print the gene-gene pair and save the ID.
                        if (!( ID in data )) {

                                print ID, ID

                                data[ID] = 1
                        }

                }
        ' $1

fi
