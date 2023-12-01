process CHECK_INPUTS {

	input:
		path gendir
		path tree
		
	output:
		path "genome_files.txt", emit : genome_files
		path "valid_tree.nwk",	 emit : tree

	"""
	chmod -R +x ${projectDir}/bin

	# 1. Check that there is a single non-empty valid FASTA and GFF3 file per genome :

	# The list of non-empty files from 'gendir' :
	gendirlist=\$(find ${gendir}/. -maxdepth 1 -not -empty -ls | awk '{print \$NF}')

	echo "\${gendirlist}" | grep "\\.fas\$\\|\\.fna\$\\|\\.fasta\$\\|\\.fa\$" | sed -E "s~.*/(.*)\\..*~\\1~" | sort | uniq -c > nb_fasta.txt
	echo "\${gendirlist}" | grep "\\.gff\$\\|\\.gff3\$" | sed -E "s~.*/(.*)\\..*~\\1~" | sort | uniq -c > nb_gff.txt
	
	awk '
		# Record every genome 
		{ gn[\$2] = 1 }

		# For any line, if the nb of occurences if not one, record the genome in "pb"
		\$1 != 1 { pb[\$2]=1 }

		# Fasta list
		FNR==NR { fas[\$2] = 1 }
		# GFF list
		FNR!=NR { gff[\$2] = 1 }

		END {
			for(file in fas){ if (!(file in gff)) { pb[file] = 1 } }
			for(file in gff){ if (!(file in fas)) { pb[file] = 1 } }

			for(file in gn){ 
				if (file in pb) { print file >> "problematic_genomes.txt" }
				else { print file >> "valid_genomes.txt" } 
			}
		}
		' nb_fasta.txt nb_gff.txt
	
	if [ -s problematic_genomes.txt ]
	then
		echo "WARNING : Each genome must be provided with one genomic FASTA file ('.fna','.fa','.fas','.fasta') and one GFF3 annotation file ('.gff3','.gff')"
		echo "ERROR : some genomes seem to not have a single non-empty valid FASTA and GFF3 file (see 'problematic_genomes.txt', 'nb_fasta.txt', 'nb_gff.txt') :"
		cat problematic_genomes.txt
		exit 1
	fi

	# The list of the genomes to process is generated from the list of FASTA files in gendir.
	echo "Genomes to procces (based on FASTA files in $gendir):"
	cat valid_genomes.txt
	echo ""


	# 2. Check that GFF files are GFF3 compliant and have "CDS", "gene", and "mRNA" features.
	cat valid_genomes.txt | while read genome
	do
		gff=$gendir/\$( ls $gendir | grep \$genome | grep "\\.gff\$\\|\\.gff3\$" )
		check_gff.sh \$gff >> gff_checking.txt
	done
	# If one or several GFF files are not correct :
	if grep -q ERROR gff_checking.txt
	then
		echo "ERROR : some GFF files are not formated as expected (see 'gff_checking.txt') :"
		cat gff_checking.txt
		exit 1
	fi


	# 3. Check the provided tree :
	if [ $tree == "dummy" ]
	then
		# If '--tree' is not provided, make a dummy empty file.
		touch valid_tree.nwk
	else
		if [ -s $tree ]
		then 
			# Check that every genome in the list is in the phylogenetic tree. If not exit.
			# Also do the reverse and exit with respect to the prune_extra_taxa argument.
			# If set to 'True', a corrected tree will be generated, with only taxa also contained in the list. Otherwise, exit. 
			# (it can be convenient to download a tree with more genome than those of immediat interest).
			check_list_n_tree.py -names valid_genomes.txt -tree $tree -prune_extra_taxa True -out valid_tree.nwk
		else
			echo "ERROR : $tree DOES NOT EXISTS or IS EMPTY."
			exit 1
		fi
	fi


	# 4. Built a channel for the rest of the pipeline
	cat valid_genomes.txt | while read genome
	do
		fasta=$gendir/\$( ls $gendir | grep \$genome | grep "\\.fas\$\\|\\.fna\$\\|\\.fasta\$\\|\\.fa\$" )
		gff=$gendir/\$( ls $gendir | grep \$genome | grep "\\.gff\$\\|\\.gff3\$" )
		echo "\${PWD}/\${fasta}__,__\${PWD}/\${gff}" >> genome_files.txt
	done

	"""


}


process EXTRACT_CDS {

	input:
		tuple val(name), path(fasta), path(gff)
	output:
		tuple val(name), path(fasta), path(gff), path("${fasta}.fai"), path("${name}_CDS.fna"), path("${name}_CDS.faa")

	"""
	chmod -R +x ${projectDir}/bin

	# Use GffRead.

	# Need a fai index
	samtools faidx $fasta

	# Keep GFF sequences that are also in the FASTA file.
	awk 'BEGIN{FS=OFS="\t"} {print "^"\$1,""} END {print "^#"}' ${fasta}.fai | grep -f - $gff > gff_filterA

	# Also remove mRNA with undefined strand and features with abnormal end
	awk -F"\t" '
			FNR==NR {max[\$1]=\$2}
			FNR!= NR && ( /^#/ || (\$7 !~ /?/ && \$5 <= max[\$1]) )
			' ${fasta}.fai gff_filterA > gff_filterB

	# Extract the genomic CDS.
	# -V discard any mRNAs with CDS having in-frame stop codons
	# -x : write a fasta file with spliced CDS for each GFF transcript
	gffread -V -g $fasta -x ${name}_CDS.fna gff_filterB


	# Linearize the .fna file (credits to https://gist.github.com/lindenb/2c0d4e11fd8a96d4c345).
	linearizefasta.sh ${name}_CDS.fna > ${name}_CDS.fna_linear
	mv ${name}_CDS.fna_linear ${name}_CDS.fna
	

	# Get the translated CDS FASTA.
	translate.sh ${name}_CDS.fna

	"""
}


process TAXDUMP {

	debug true

	input:
		path taxdump
		
	output:
		env valid_taxdump
		
	"""
	# If '--taxdump' is not provided,
	if [ $taxdump == "dummy" ]
	then
		if [ ! -d ${workDir}/taxdump ]
		then
			# Download it
			echo "Downloading ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
			wget -N ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
			mkdir ${workDir}/taxdump && tar zxf new_taxdump.tar.gz -C ${workDir}/taxdump
		fi
		valid_taxdump="${workDir}/taxdump"
		echo "The taxdump directory is in \${valid_taxdump}. You can use '--taxdump \${valid_taxdump}'."

	else

		# Just keep taxdump as it is.
		echo "The user provided a taxdump directory. Using it."
		valid_taxdump=$taxdump

	fi
	"""
}


process GENERA {

	tag 'phylostratigraphy'

	publishDir "${params.outdir}/genera_results"

	cpus params.max_cpus

	input:
		val focal_taxid
		path focal_CDS
		path neighbors_CDS
		path db
		path taxdump
		
	output:
		path "${focal_taxid}_gene_ages.tsv"
		
	"""
	# -t taxID
	# -q query fasta
	# -n number of CPU to use
	# -b path to the nr database
	# -r OR, the first time you use genEra : -d /path/to/taxdump/ \
	# -x Temp dir (potentially hundreds of Go needed).
	genEra \
	-t $focal_taxid \
	-q $focal_CDS \
	-a $neighbors_CDS \
	-n ${task.cpus} \
	-b $db/nr \
	-d $taxdump
	"""
}


process GENERA_FILTER {

        publishDir "${params.outdir}/genera_results"

        input:
            path taxdump
            val taxid
            val TRG_node
            val TRG_rank
            path mapping
            path genera_out

        output:
            path 'TRG_CDS.txt'

        """
        genEra_filter.py '$taxdump' '$taxid' '$TRG_node' '$TRG_rank' '$mapping' '$genera_out'
        """
}


process FIND_TRG {

	debug true

	input:
		path(TRG_list)
		tuple val(genome_name), path(gfasta), path(gff), path(fai), path(CDS_fna), path(CDS_faa)
	output:
		path("TRG.txt")
		
	"""
	# Get the CDS list in a text file (already unique with gffread).
	grep ">" $CDS_fna | sed "s/>//" > CDS.txt


	if [ -s $TRG_list ]
	then
		# Only keep TRG that are present in the CDS exctracted by gffread.
		grep -x -f CDS.txt $TRG_list | sort | uniq > TRG.txt || true

		if [ -s TRG.txt ]
		then
			echo "Processing \$(cat TRG.txt | wc -l) TRGs..."

			grep -v -x -f CDS.txt $TRG_list | sort | uniq > TRG_not_exctracted_by_gffread.txt || true

			if [ -s TRG_not_exctracted_by_gffread.txt ]
			then
				echo "WARNING : some provided TRGs were not exctracted by gffread (see 'TRG_not_exctracted_by_gffread.txt')."
				cat TRG_not_exctracted_by_gffread.txt
			fi
		else
			echo "ERROR : none of the provided TRGs was exctracted by gffread!\nThe TRGs list must contain CDS names (in the GFF3 file : 'mRNA' features' ID)."
			exit 1
		fi
	else
		echo "ERROR : ${TRG_list} does not exist or is empty."
		exit 1
	fi
	"""
}


process TRG_FNA {

	input:
		path CDS_fna
		path TRGs
		
	output:
		path "TRG.fna"
		
	"""
	chmod -R +x ${projectDir}/bin

	faSomeRecords $CDS_fna $TRGs TRG.fna
	"""
}


process MULTIELONGATE_FOCAL_TRG {

	input:
		path TRGs
		tuple val(genome), path(gfasta), path(gff), path(fai), path(CDS_fna, stageAs: "CDS.fna"), path(CDS_faa, stageAs: "CDS.faa")
		
	output:
		path "TRG_multielongated.faa", emit : TRG_multielongated_faa
		
	"""
	chmod -R +x ${projectDir}/bin

	# Generate a FASTA with the CDS of the focal genome, elongated by 100 nucleotides up and downstream :
	
	# 99_nucl_translated_in_frame_0____translated_CDS____99_nucl_translated_in_frame_0
	# 99_nucl_translated_in_frame_0____translated_CDS____99_nucl_translated_in_frame_1
	# 99_nucl_translated_in_frame_0____translated_CDS____99_nucl_translated_in_frame_2
	# 99_nucl_translated_in_frame_1____translated_CDS____99_nucl_translated_in_frame_0
	# 99_nucl_translated_in_frame_1____translated_CDS____99_nucl_translated_in_frame_1
	# 99_nucl_translated_in_frame_1____translated_CDS____99_nucl_translated_in_frame_2
	# 99_nucl_translated_in_frame_2____translated_CDS____99_nucl_translated_in_frame_0
	# 99_nucl_translated_in_frame_2____translated_CDS____99_nucl_translated_in_frame_1
	# 99_nucl_translated_in_frame_2____translated_CDS____99_nucl_translated_in_frame_2
	
	extend_cds.py \
	--gfasta $gfasta \
	--fai $fai \
	--gff $gff \
	--cds $TRGs \
	--size 99 \
	--mode multi \
	--out TRG_multielongated.faa
	"""
}


process ELONGATE_CDS {

	input:
		tuple val(name), path(gfasta), path(gff), path(fai), path(CDS_fna, stageAs: "CDS.fna"), path(CDS_faa, stageAs: "CDS.faa")

		
	output:
		tuple val(name), path(gfasta), path("CDS_elongated.faa")
		
	"""
	chmod -R +x ${projectDir}/bin

	## In order to search for homologous CDS in the genome, get a FASTA with the elongated CDS of the genome (100 nucl upstream and downstream).
	# In case of multi matchs, the one with the best evalue is more likely to be the right one with this elongated version (bring genomic context).
	
	# Generate the elongated translated CDS FASTA.
	extend_cds.py \
	--gfasta $gfasta \
	--fai $fai \
	--gff $gff \
	--cds $CDS_fna \
	--size 99 \
	--mode simple \
	--out CDS_elongated.faa
	"""
}


process BLAST {

	label 'parallel_job'

	publishDir "${params.outdir}/blast_out", mode: 'copy',  pattern: "*.out"

	cpus   params.max_cpus

	input:
		path(TRG_multielongated_faa)
		tuple val(genome_name), path(gfasta), path(CDS_elongated_faa)
	output:
		tuple val(genome_name), path("*_CDS_elongated.out"), path("*_genome.out")
		
	"""

	# Create the required BLAST databases :
	makeblastdb -dbtype prot -in $CDS_elongated_faa -out ${CDS_elongated_faa}.db
	
	# Also create a nucleotide db with the entiere genome, to look for any matchs (only considered if no CDS matchs are found).
	makeblastdb -dbtype nucl -in $gfasta -out ${gfasta}.db
		

	# Perform a ### BLASTp ### of the multielongated translated TRG FASTA of the focal genome against the elongated CDS FASTA of the current genome.
	# The -num_threads option should be adapted.
	# Header = "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs sframe"
	blastp -query $TRG_multielongated_faa -db ${CDS_elongated_faa}.db -num_threads ${task.cpus} \
	-outfmt "7 std qlen qcovs sframe" \
	-evalue 0.001 > TRG_multielongated_blastp_${genome_name}_CDS_elongated.out
	

	# Similarly, perform a ### tBLASTn ### of the multielongated translated TRG FASTA of the focal genome against the genomic nucleotide FASTA of the current genome.
	# Header = " qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp sframe"
	tblastn -query $TRG_multielongated_faa -db ${gfasta}.db -num_threads ${task.cpus} \
			-outfmt "7 std qlen qcovhsp sframe" \
	-evalue 0.001 > TRG_multielongated_tblastn_${genome_name}_genome.out
	
	"""
}


process TRG_LIST_BEFORE_STRATEGY {

	input:
		tuple val(genome_name), path(TRG_multielongated_blastp_CDS_elongated_out), path(TRG_multielongated_tblastn_genome_out)
		path focal_mRNA_to_gene
		
	output:
		path 'TRG_before_strategy.tsv'
		
	"""
	# Get the list of TRGs before applying the selection strategy.
	if grep -q "# Query:" $TRG_multielongated_blastp_CDS_elongated_out
	then
		# BLAST format 7
		awk '
			/# Query:/ {
				CDS=gensub(/# Query:\s(.*)/,"\\\\1","g",\$0)
				# If the case where the CDS is elongated
				short=gensub(/_elongated.*/,"","g",CDS)
				if(!(short in data)){
					print short
					data[short]=1
				}
			}
		' $TRG_multielongated_blastp_CDS_elongated_out > TRG_CDS_before_strategy.txt
	else
		# BLAST format 6 (diamond includes every query including those without any match) ## NOT IMPLEMENTED ##
		awk -F"\t" '
			\$0 !~ /^#/ {
				short=gensub(/_elongated.*/,"","g",\$1)
				if(!(short in data)){
					print short
					data[short]=1
				}
			}
		' $TRG_multielongated_blastp_CDS_elongated_out > TRG_CDS_before_strategy.txt

	fi

	awk '
		BEGIN {FS=OFS="\t"}

		# First file : focal_mRNA_to_gene
		FNR==NR { genes[\$1]=\$2 }

		# Second file : TRGs_before_strategy
		FNR!=NR {

			if(\$1 in genes){
				print \$0, genes[\$1]
			} else {
				if(FNR != 1) { print "WARNING : NO gene associated with "\$1"!" }
			}
		}
	' $focal_mRNA_to_gene TRG_CDS_before_strategy.txt > TRG_before_strategy.tsv
	"""
}


process BLAST_FILTER {

	input:
		val focal
		tuple val(genome_name), path(TRG_multielongated_blastp_CDS_elongated_out), path(TRG_multielongated_tblastn_genome_out)
		
	output:
		tuple val(genome_name), path('TRG_blast_*_best_hits.tsv')
		
	"""
	chmod -R +x ${projectDir}/bin

	# Get the list of homologs (CDS and whole_genome) for each genome.
	
	bash tableau_frameshifts_pretreatment.sh --in $TRG_multielongated_tblastn_genome_out
	
	bash tableau.sh --query_genome $focal --subject_genome $genome_name --in $TRG_multielongated_blastp_CDS_elongated_out 		--type CDS
	bash tableau.sh --query_genome $focal --subject_genome $genome_name --in ${TRG_multielongated_tblastn_genome_out}.processed 	--type genome
	
	"""
}


process DUMMY_DISTANCES {

	input:
		path(names)
	output:
		path "dummy_distances.tsv"

	"""
	awk '
		BEGIN {
			FS=OFS="\t"
			print "Taxon","Time of divergence (MYA)","Root distance"	
			}
		{ print \$0,"NA","NA" }

	' $names > dummy_distances.tsv
	"""
}


process TREE_DISTANCES {

	input:
		val focal
		path tree
	output:
		path "tree_distances.tsv"

	"""
	chmod -R +x ${projectDir}/bin

	tree_distances.py -tree $tree -focal $focal -out tree_distances.tsv 
	"""
}


process TRG_TABLE {

	input:
		path tree_distances
		path best_hits
		
	output:
		path 'TRG_table.tsv'
				
	"""
	# Generate a TRG table that will be easily accessible to the user.
	# It will serve as raw material to select good candidates.
	awk '
	
		BEGIN { 
			FS=OFS="\t"
			print "Focal","TRG","genome","sequence","type","time of divergence (MY)","root distance","isSynt"	
		}

		# First file (i.e. tree_distances)
		# Save time of divergence in TD and root distance in RD.
		FNR==NR { TD[\$1]=\$2 ; RD[\$1]=\$3	}

		# Second file (i.e. TRG_best_hits )
		# Print the line and add the time of divergence, the root distance, and "NA"
		# This last column will be filled with the CHECK_SYNTENY ouputs if appropriate.
		FNR!=NR { print \$0, TD[\$3], RD[\$3], "NA" }

	' $tree_distances $best_hits > TRG_table.tsv
	"""
}


process CHECK_SYNTENY_INPUTS {

	input:
		val (focal)
		tuple val(genome), path(best_hits)
		
	output:
		tuple val(genome), path("${focal}_vs_${genome}_synteny_pairs.tsv")
		
	"""
	touch ${focal}_vs_${genome}_synteny_pairs.tsv
	awk '
		BEGIN {FS=OFS="\t"}
		
		\$5 == "CDS" { data[\$2\$3] = 1 }
		\$5== "genome" && (!(\$2\$3 in data)) { 
			print gensub(/_elongated.*/,"","g",\$2),gensub(/(.*)_([0-9]+)_([0-9]+)\$/,"\\\\1\t\\\\2\t\\\\3","g",\$4) }' $best_hits > ${focal}_vs_${genome}_synteny_pairs.tsv 
	"""
}


process CHECK_SYNTENY {

	publishDir "${params.outdir}/synteny"

	// max_proc_mem = '20.GB'
	// memory { "${MemoryUnit.of(params.max_memory).toMega()}mb" < "${MemoryUnit.of(max_proc_mem).toMega()}mb" ? params.max_memory : max_proc_mem }
	
	input:
		val anch_nb
		tuple val(focal), path(focal_gff, stageAs:'focal_gff.gff')
		tuple val(genome), path(genome_gff), path(orthologs), path(input_pairs)
		
	output:
		tuple path("${focal}_vs_${genome}_synteny.out"), path("${focal}_vs_${genome}_synteny.tsv")
		
	"""
	chmod -R +x ${projectDir}/bin

	python --version
	touch ${focal}_vs_${genome}_synteny.out ${focal}_vs_${genome}_synteny.tsv
	check_synteny.py \
	--gffA $focal_gff \
	--gffB $genome_gff \
	--ortho $orthologs \
	--list $input_pairs \
	--flankA $anch_nb \
	--flankB $anch_nb \
	--up_min 1 \
	--ov_min 0 \
	--do_min 1 \
	--insertion 0 \
	--ins_max -1 \
	--del_max -1 \
	--rea_max -1 \
	--erase \
	--out ${focal}_vs_${genome}_synteny
	"""
}


process SYNTENY_TO_TABLE {

	publishDir "${params.outdir}"
	
	input:
		path TRG_table
		path syntenys

	output:
		path "TRG_table.tsv"
		
	"""
	awk -v TRG_table="${TRG_table}" '

		BEGIN {FS=OFS="\t"}

		# In the synteny files, except for the header,
		# store data[query_vs_subject] = isSynt
		FNR!= 1 && FILENAME != TRG_table {

			# synteny files are named "focalname_vs_neighborname_synteny.tsv"
			neighbor = gensub(/.*_vs_(.*)_synteny.tsv/,"\\\\1","g",FILENAME)
			query  = \$1
			target = \$2
			isSynt = \$9

			data[query"_vs_"neighbor":"target]= isSynt

		}

		# In the TRG table,
		# print the header
		FILENAME == TRG_table && FNR == 1

		# then
		FILENAME == TRG_table && FNR != 1 {

			query = gensub(/_elongated.*/,"","g",\$2)
			neighbor = \$3
			target = \$4

			# If the query CDS has been tested against the given target region in the given neighbor,
			if(query"_vs_"neighbor":"target in data){
				# Fill the "isSynt" column with True or False (the checking result).
				\$8 = data[query"_vs_"neighbor":"target]
			}

			print \$0

		}
	' $syntenys $TRG_table > ${TRG_table}.tmp

	mv ${TRG_table}.tmp ${TRG_table}
	"""
}


process FILTER_TABLE_WITH_STRATEGY {

	input:
		val focal
		val synteny
		val strategy
		path TRG_table
		
	output:
		path "TRGs_selected_before_isoform_control.txt"
		
	"""
	filter_table_with_strategy.py --table $TRG_table --strategy $strategy --synteny $synteny --out TRGs_selected_before_isoform_control.txt
	"""
}


process FILTER_ISOFORMS {

	publishDir "${params.outdir}"

	input:
		path TRGs_selected_before_strategy
		path TRGs_selected_after_strategy
		
	output:
		path "denovogenes.txt"
		
	"""
	touch denovogenes.txt

	# Compare the list of TRGs before and after applying the slection strategy.
	# Remove the TRGs with excluded isoforms.
	awk '
		BEGIN {FS=OFS="\t"}

		FNR==1 { FNUM++ }

		# First file : TRGs_selected_after_strategy
		FNUM==1 { final[\$0]=1 }
		
		# Second file : TRGs_selected_before_strategy
		FNUM==2 && (!(\$1 in final)){ excluded[\$2]=1 }
		
		# Third file : TRGs_selected_before_strategy
		FNUM==3 {
			if(\$2 in excluded){
				print \$1 >> "TRG_excluded.txt"
			}
			else {
				print \$0 >> "denovogenes.txt"
			}
		}

	' $TRGs_selected_after_strategy $TRGs_selected_before_strategy $TRGs_selected_before_strategy
	"""
}


process WARN_MISSING {
	debug true
	input:
	val x
	val y

	when:
    y == 'EMPTY'

	"""
	echo $x
	exit 1
	"""
}
