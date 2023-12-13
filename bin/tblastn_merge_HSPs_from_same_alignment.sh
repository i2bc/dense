#!/bin/bash

# The idea behind this script is to identify alignments made of several HSP due to frameshifts in the subject.
# When BLAST considers two or more HSP to be part of a same alignment, it recomputes a new evalue with sum-statistics.
# all the different HSP of the same allignment thus share the same evalue but keep their own bitscore.

# The tblastn command has been configured with -outfmt "7 std qlen qcovhsp sframe"

# Get bash named parameters (adapted from https://www.brianchildress.co/named-parameters-in-bash/).
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
  shift
done




awk '

        BEGIN{ FS=OFS="\t" }


        # For non-comment lines (the evalue is assumed to have been set in the blast command),
        $1 !~ /#/ {

                qlen=$13

                # Get the current subject strand
                if( $15 > 0){   curr_sstrand= "+" }
                else {          curr_sstrand= "-" }

                # If the query and subject are the same as for the previous HSP (line), and the evalue is non-null and the same as the previous HSP, and with the same strand but with a different bitscore,
                # then is it assumed that the current HSP belongs to the same alignment as the previous one.
                if( $1 == query && $2 == subject && $11 != 0 && $11 == evalue && curr_sstrand == sstrand && $12 != bitscore ){

                        new_qs=$7
                        new_qe=$8
                        new_ss=$8
                        new_se=$9

                        # Get the earliest query-start and the lattest query end among the current and the previous HSP.
                        if( new_qs > qs) { new_qs=qs }
                        if( new_qe < qe) { new_qe=qe }

                        # With repesct to the subject strand, also elongate the current HSP if appropriate.
                        # Positive strand
                        if( sstrand ~ /\+/ ){
                                if( new_ss > ss) { new_ss=ss }
                                if( new_se < se) { new_se=se }
                                HSP_slen= 1 + new_se - new_ss
                        }
                        # Negative strand
                        if( sstrand ~ /\-/ ){
                                if( new_ss < ss) { new_ss=ss }
                                if( new_se > se) { new_se=se }
                                HSP_slen= 1 + new_ss - new_se
                        }

                        # If the potential elongated HSP is no more than 3 times (3 times fors AA, 9 for nucl) greater than the query, apply the modifications (e.i. perfom the elongation/merging).
                        if( HSP_slen <= 9*qlen ){
                                $7=new_qs
                                $8=new_qe
                                $9=new_ss
                                $10=new_se
                        }
                }

                # Save the current HSP features for next comparison.
                query=$1
                subject=$2
                qs=$7
                qe=$8
                ss=$9
                se=$10
                evalue=$11
                bitscore=$12
                sstrand= curr_sstrand

                # Recompute the qcovhsp field
                HSP_qlen= 1 + qe - qs
                $14= 100 * HSP_qlen / qlen
                qcov=$14

                # Rename the subject field to include the HSP"s subject start and end positions.
                $2=$2"_"ss"_"se

                # Print the (possibly edited) current HSP.
                print $0

        }

' $in > ${in}.processed
