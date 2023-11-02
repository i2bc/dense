process BLAST {

	label 'parallel_job'

	publishDir "${params.outdir}/blast_out", mode: 'copy',  pattern: "*.out"

	cpus params.max_cpus
	
	input:
		tuple val(focal_name), path(focal_CDS_faa, stageAs: 'focal_CDS.faa')
		tuple val(neighbor_name), path(neighbor_CDS_faa)

	output :
		tuple val(focal_name), val(neighbor_name), path("${focal_name}_CDS_BLASTp_${neighbor_name}_CDS*.out"), path("${neighbor_name}_CDS_BLASTp_${focal_name}_CDS*.out")

	"""
	chmod -R +x ${projectDir}/bin

	echo ${task.cpus}
	timestamp=\$(date -d "today" +"%Y%m%d_%H%M")

	# Create the required BLAST databases :
	makeblastdb -dbtype prot -in $neighbor_CDS_faa -out $neighbor_CDS_faa


	# Perform a ### BLASTp ### of the translated CDS FASTA of the focal genome against the translated CDS FASTA of the current genome.
	# The -num_threads option should be adapted.
	# blastp header = "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs"
	blastp -query $focal_CDS_faa -db $neighbor_CDS_faa -num_threads ${task.cpus} \
	-outfmt "7 std qlen qcovs" \
	-evalue 0.001 > ${focal_name}_CDS_BLASTp_${neighbor_name}_CDS_\${timestamp}.out	
	

	# Do the reciprocal (except when focal and neighbor are the same).
	if [ ${focal_name} != ${neighbor_name} ]
	then
		makeblastdb -dbtype prot -in $focal_CDS_faa -out $focal_CDS_faa
		
		# blastp header = "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs"
		blastp -query $neighbor_CDS_faa -db $focal_CDS_faa -num_threads ${task.cpus} \
		-outfmt "7 std qlen qcovs" \
		-evalue 0.001 > ${neighbor_name}_CDS_BLASTp_${focal_name}_CDS_\${timestamp}.out	
	fi 
	"""
}


process DIAMOND_BLAST {

	label 'parallel_job'

	publishDir "${params.outdir}/diamondblast_out", mode: 'copy',  pattern: "*.out"
	
	cpus params.max_cpus

	input:
		tuple val(focal_name), path(focal_CDS_faa, stageAs: 'focal_CDS.faa')
		tuple val(neighbor_name), path(neighbor_CDS_faa)
		val sensitivity

	output :
		tuple val(focal_name), val(neighbor_name), path("${focal_name}_CDS_BLASTp_${neighbor_name}_CDS*.out"), path("${neighbor_name}_CDS_BLASTp_${focal_name}_CDS*.out")

	"""
	chmod -R +x ${projectDir}/bin

	echo "cpus = ${task.cpus}"
	timestamp=\$(date -d "today" +"%Y%m%d_%H%M")

	# Create the required databases :
	diamond makedb --in $neighbor_CDS_faa -d $neighbor_CDS_faa --threads ${task.cpus}


	# Perform a ### BLASTp ### of the translated CDS FASTA of the focal genome against the translated CDS FASTA of the current genome.
	# The -num_threads option should be adapted.
	# blastp header = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
	diamond blastp  --query $focal_CDS_faa --db $neighbor_CDS_faa --threads ${task.cpus} \
	-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
	--evalue 0.001 --unal 1 -k 0 --max-hsps 0 \
	-o ${focal_name}_CDS_BLASTp_${neighbor_name}_CDS_\${timestamp}.out \
	$sensitivity # optional sensitivity setting

	add_qcovs.sh ${focal_name}_CDS_BLASTp_${neighbor_name}_CDS_\${timestamp}.out > ${focal_name}_CDS_BLASTp_${neighbor_name}_CDS_\${timestamp}.out.tmp
	mv ${focal_name}_CDS_BLASTp_${neighbor_name}_CDS_\${timestamp}.out.tmp ${focal_name}_CDS_BLASTp_${neighbor_name}_CDS_\${timestamp}.out


	# Do the reciprocal (except when focal and neighbor are the same).
	if [ ${focal_name} != ${neighbor_name} ]
	then
		diamond makedb --in $focal_CDS_faa -d $focal_CDS_faa --threads ${task.cpus}
		
		# blastp header = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
		diamond blastp  --query $neighbor_CDS_faa --db $focal_CDS_faa --threads ${task.cpus} \
		-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
		--evalue 0.001 --unal 1 -k 0 --max-hsps 0 \
		-o ${neighbor_name}_CDS_BLASTp_${focal_name}_CDS_\${timestamp}.out \
		$sensitivity # optional sensitivity setting

		add_qcovs.sh ${neighbor_name}_CDS_BLASTp_${focal_name}_CDS_\${timestamp}.out > ${neighbor_name}_CDS_BLASTp_${focal_name}_CDS_\${timestamp}.out.tmp
		mv ${neighbor_name}_CDS_BLASTp_${focal_name}_CDS_\${timestamp}.out.tmp ${neighbor_name}_CDS_BLASTp_${focal_name}_CDS_\${timestamp}.out
	fi 
	"""
}


process BEST_HITS {
	
	input:
		tuple val(focal_name), val(neighbor_name), path( focal_vs_neighbor_out, stageAs: 'focal_neighbor.out' ), path( neighbor_vs_focal_out )
		
	output:
		tuple val(neighbor_name), path('focal_neighbor_best_hits.tsv'), path("${neighbor_vs_focal_out.getBaseName()}_best_hits.tsv"), emit : mainout
		val focal_name, emit : focal_name

		
	"""
	for file in $focal_vs_neighbor_out $neighbor_vs_focal_out
	do
		awk '
		BEGIN {FS=OFS="\t"}
		
		# For HSP (non-commment) lines,
		# with a query coverage by subject >= 70%, and an unprinted query,
		# Save the query and print the query-subject pair.
		\$0 !~ /#/ && \$14 >= 70 && (!(\$1 in data )) { data[\$1]=1 ; print \$1,\$2 }
		
		' \$file > \$(basename \$file .out)_best_hits.tsv
	done
	"""
}


process MRNA_TO_GENE {
	
	input:
		tuple val(name), path(gff)
		
	output:
		tuple val(name), path("${name}_mRNA_to_gene.tsv")
		
	"""
	chmod -R +x ${projectDir}/bin

	get_mRNA_parent.sh $gff > ${name}_mRNA_to_gene.tsv      
	"""
}


process RECIPROCAL_BEST_HITS {

	publishDir "${params.outdir}/orthologs", mode: 'copy'
	
	input:
		tuple val(focal_name), path(focal_mRNA_to_gene, stageAs : 'focal_mRNA_to_gene.tsv')
		tuple val(neighbor_name), path(focal_neighbor_best_hits), path(neighbor_focal_best_hits), path(neighbor_mRNA_to_gene)
		
	output:
		tuple val(neighbor_name), path("${focal_name}_${neighbor_name}_orthologs.tsv")
		
	"""
	touch ${focal_name}_${neighbor_name}_orthologs.tsv

	# Turn A_CDS_bh_B_CDS and B_CDS_bh_A_CDS to A_gene_bh_B_gene and B_gene_bh_A_gene_REV (the rev is to be able to search for exacte same lines)
	awk '
		BEGIN {FS=OFS="\t"}

		FNR==1 {FNUM++}

		FNUM==1 && NF==2 { A_gene[\$1]=\$2 }

		FNUM==2 && NF==2 { B_gene[\$1]=\$2 }

		FNUM==3 && NF==2 && (\$1 in A_gene) && (\$2 in B_gene) { print A_gene[\$1],B_gene[\$2] > "A_gene_vs_B_gene"     }

		FNUM==4 && NF==2 && (\$2 in A_gene) && (\$1 in B_gene) { print A_gene[\$2],B_gene[\$1] > "B_gene_vs_A_gene_REV" }

	' $focal_mRNA_to_gene $neighbor_mRNA_to_gene $focal_neighbor_best_hits $neighbor_focal_best_hits

	# Search for exact same lines, and get unique lines
	grep -Fxf A_gene_vs_B_gene B_gene_vs_A_gene_REV | \
	awk -F"\t" '(!(\$0 in data)) { print \$0 ; data[\$0]=1 }' > ${focal_name}_${neighbor_name}_orthologs.tsv
	"""
}
