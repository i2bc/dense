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

                        # With repesct to the subject strand, aslo elongate the current HSP if appropriate.
                        if( sstrand ~ /\+/ ){
                                if( new_ss > ss) { new_ss=ss }
                                if( new_se < se) { new_se=se }

                                # If the potential elongated HSP is no more than 3 times (3 times fors AA, 9 for nucl) greater than the query, apply the modifications.
                                if( (1 + new_se - new_ss) <= 9*$13 ){
                                        $7=new_qs
                                        $8=new_qe
                                        $9=new_ss
                                        $10=new_se
                                }
                        }
                        if( sstrand ~ /\-/ ){
                                if( new_ss < ss) { new_ss=ss }
                                if( new_se > se) { new_se=se }
                                if( (1 + new_ss - new_se) <= 9*$13 ){
                                        $7=new_qs
                                        $8=new_qe
                                        $9=new_ss
                                        $10=new_se
                                }
                        }

                }

                # Recompute the qcovhsp field
                $14= 100 * (1 + $8 - $7) / $13

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
                qcov=$14


                $2=$2"_"ss"_"se

                # Print the (possibly edited) current HSP.
                print $0

        }

' $in > ${in}.processed
