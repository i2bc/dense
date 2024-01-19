awk -F"\t" '

BEGINFILE { printf FILENAME" : " }

$0 !~ /^#/ && $3 ~ /^CDS$/ && $9 ~ /Parent=/ { has_CDS = 1 }

$0 !~ /^#/ && $3 ~ /^mRNA$/ && $9 ~ /ID=/ && $9 ~ /Parent=/ { has_mRNA = 1 }

$0 !~ /^#/ && $3 ~ /^gene$/ && $9 ~ /ID=/ { has_gene = 1 }

has_mRNA && has_gene { valid = 1 ; exit 0 }

END {
	if( valid ) { print "valid format" }
	else {
		print "ERROR - invalid format :"
		print FILENAME >> "gff_to_correct.txt"
		if (!(has_CDS)) { print "The GFF file must have \"CDS\" features with a \"Parent\" attribute (9th column)." }
		if (!(has_mRNA)){ print "The GFF file must have \"mRNA\" features with both \"ID\" and \"Parent\" attributes (9th column)." }
		if (!(has_gene)){ print "The GFF file must have \"gene\" features with an \"ID\" attribute (9th column)." }
	}
}

' $1
