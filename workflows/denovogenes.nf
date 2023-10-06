/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WELCOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Here is an attempt of ASCII art
log.info "Project : $workflow.projectDir"
log.info "      _                                                           "
log.info "     | |                                                          "
log.info "   __| | ___   _ __   _____   _____     __ _  ___ _ __   ___  ___ "
log.info "  / _` |/ _ \\ | '_ \\ / _ \\ \\ / / _ \\   / _` |/ _ \\ '_ \\ / _ \\/ __|"
log.info " | (_| |  __/ | | | | (_) \\ V / (_) | | (_| |  __/ | | |  __/\\__ \\"
log.info "  \\__,_|\\___| |_| |_|\\___/ \\_/ \\___/   \\__, |\\___|_| |_|\\___||___/"
log.info "                                        __/ |                     "
log.info "                                       |___/                      "
log.info ""
log.info "Welcome to the de novo genes Nextflow pipeline."
log.info "\n"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS MANAGMENT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run denovogenes.nf --gendir <DIR WITH GFF AND FASTA> --focal <FOCAL_GENOME_NAME> --outdir <OUTDIR>")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


// Parameters compatibility
stop = false

if (params.strategy == 1 && !params.tree) {
	log.info "A phylogenetic tree must be supplied with strategy 1 (see '--strategy' and '--tree')."
	stop = true
}

if (params.strategy != 1 && params.tree) {
	log.info "WARNING : the provided phylogenetic tree will not be used with strategy ${params.strategy} (see '--strategy' and '--tree')."
}

if (params.trg_node && params.trg_rank) {
	log.info "'--trg_node' and '--trg_rank' are not compatible'. Select only one of them."
	stop = true
}

if (params.trgs) {
	for (key in ['genera_out', 'genera_db', 'taxids', 'trg_node', 'trg_rank']) {
		if(params[key]){
			log.info "'--trgs' and '--${key}' are not compatible."
			stop = true
		}
	}
}

if (params.genera_out) {
	if(!params.taxids){
		log.info "'--taxids' is required with '--genera_out'."
		stop = true
	}
	if (!params.trg_node && !params.trg_rank) {
		log.info "'--trg_rank' OR '--trg_node' is required with '--genera_out'."
		stop = true
	}
}

if ( !params.taxids && params.genera_db ) {
	log.info "'--taxids' is required with '--genera_db."
	stop = true
}

def synt_strats = [1, 2]
if (params.synteny && ! synt_strats.contains(params.strategy)) {
	log.info "There is no synteny check to perform with strategy ${params.strategy} since there no homolog detection (see '--strategy' and '--synteny')."
	stop = true
}

for (key in ['orthodir', 'blastdir']) {
	if (params[key] && ! synt_strats.contains(params.strategy)) {
		log.info "There is no synteny check to perform with strategy ${params.strategy} since there no homolog detection (see '--strategy' and '--${key}')."
		stop = true
	}
	if (params[key] && ! params.synteny) {
		log.info "'--${key}' needs '--synteny' to be set to 'true'."
		stop = true
	}
}

if (stop) {
	exit 0
}


/// tree
if(params.tree){
	tree = file(params.tree)
} else {
	tree = file("dummy")
}

/// taxdump
if(params.taxdump){
	taxdump = file(params.taxdump)
} else {
	taxdump = file("dummy")
}

// trg_node
if(!params.trg_node){ trg_node = "null"} else { trg_node = params.trg_node }

// TRG list file
// if ( params.containsKey('trgs') && (params.trgs instanceof String) && params.trgs != "null") {
// 	if (file(params.trgs).exists()){
// 		if(file(params.trgs).isEmpty()){
// 			log.info "WARNING : '${params.trgs}' is empty (see the '--trgs' option)"
// 			System.exit(0)
// 		} else {
// 			log.info "A list of TRG was provided with '${params.trgs}'."
// 			TRG_file = params.trgs
// 		}
// 	} else {
// 		log.info "WARNING : '${params.trgs}' not found (see the '--trgs' option)"
// 		System.exit(0)
// 	}
// } else { 
// 	log.info "No list of TRGs provided."
// 	TRG_file = file("dummy")
// }
if (params.trgs) {
	if (file(params.trgs).exists()){
		if(file(params.trgs).isEmpty()){
			log.info "WARNING : '${params.trgs}' is empty (see the '--trgs' option)"
			System.exit(0)
		} else {
			log.info "A list of TRG was provided with '${params.trgs}'."
		}
	} else {
		log.info "WARNING : '${params.trgs}' not found (see the '--trgs' option)"
		System.exit(0)
	}
}

if ( params.blastdir ){ log.info "A directory with the precomputed BLAST outputs has been provided : '${params.blastdir}'" }
if ( params.orthodir ){ log.info "A directory with the precomputed ortholog pairs has been provided : '${params.orthodir}'" }
log.info ""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ORTHOLOGS                  } from '../subworkflows/local/orthologs'
include { CHECK_INPUTS               } from '../modules/local/denovogenes_modules.nf'
include { MRNA_TO_GENE               } from '../modules/local/orthologs_modules.nf'
include { EXTRACT_CDS                } from '../modules/local/denovogenes_modules.nf'
include { GENERA                     } from '../modules/local/denovogenes_modules.nf'
include { GENERA_FILTER              } from '../modules/local/denovogenes_modules.nf'
include { FIND_TRG                   } from '../modules/local/denovogenes_modules.nf'
include { MULTIELONGATE_FOCAL_TRG    } from '../modules/local/denovogenes_modules.nf'
include { ELONGATE_CDS               } from '../modules/local/denovogenes_modules.nf'
include { BLAST                      } from '../modules/local/denovogenes_modules.nf'
include { TRGS_BEFORE_STRATEGY       } from '../modules/local/denovogenes_modules.nf'
include { BLAST_FILTER               } from '../modules/local/denovogenes_modules.nf'
include { DUMMY_DISTANCES            } from '../modules/local/denovogenes_modules.nf'
include { TREE_DISTANCES             } from '../modules/local/denovogenes_modules.nf'
include { TRG_TABLE                  } from '../modules/local/denovogenes_modules.nf'
include { CHECK_SYNTENY_INPUTS       } from '../modules/local/denovogenes_modules.nf'
include { CHECK_SYNTENY              } from '../modules/local/denovogenes_modules.nf'
include { SYNTENY_TO_TABLE           } from '../modules/local/denovogenes_modules.nf'
include { FILTER_TABLE_WITH_STRATEGY } from '../modules/local/denovogenes_modules.nf'
include { FILTER_ISOFORMS            } from '../modules/local/denovogenes_modules.nf'
include { TEST                       } from '../modules/local/denovogenes_modules.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DENOVOGENES {
//workflow {

	// TEST()

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Check the inputs before the beginning.
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	CHECK_INPUTS(
				 file(params.gendir),
				 tree
				)

	// genomic FASTA and GFF3 pairs that will be processed
	genome_ch = CHECK_INPUTS.out.genome_files
		.splitText() 
		.map { tuple( it.strip().split("__,__") ) }
		.map { fasta, gff -> [ file(fasta).getSimpleName(), fasta, gff ] }

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Exctract the CDS of each genome
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	EXTRACT_CDS(genome_ch)

	// A channel with only the focal genome
	focal_ch = EXTRACT_CDS.out
	.map { name, fasta, gff, fai, CDS_fna, CDS_faa -> [ name, fasta, gff, fai, CDS_fna, CDS_faa ] }
	.filter { name, fasta, gff, fai, CDS_fna, CDS_faa -> name == params.focal }
	// .filter ( ~/\[${params.focal}, .*/ )

	// focal mRNA_to_gene mapping
	MRNA_TO_GENE(
				 focal_ch
				 .map{ name, fasta, gff, fai, CDS_fna, CDS_faa -> [ name, gff ] }
				)
	focal_mRNA_to_gene_ch = MRNA_TO_GENE.out
	.map{ name, mapping -> mapping }
	.first()

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Get the TRGs and their homologous sequences (CDS or genomic intervals)
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	// 'If a directory with the blast output of the TRG (queries) is provided.'
	if (params.trgsblastdir) { 

		// Skip the next block since the TRGs collection and the search for their hommologs are alreayd computed.
		// Mimic its output.
		// tuple val(genome_name), path("*_CDS_elongated.out"), path("*_genome.out")
		BLAST_out_ch = genome_ch
		.map { name, fasta, gff -> [ name,
									 file("${params.trgsblastdir}/TRG_multielongated_blastp_${name}_CDS_elongated.out"),
									 file("${params.trgsblastdir}/TRG_multielongated_tblastn_${name}_genome.out")
								   ] 
			 }

	} else {

		// First get the TRG list

		// If there is a predefined list of TRG
		if (params.trgs) {

			// Simply collect the list (with a quick check-up)
			FIND_TRG( 
					 file(params.trgs),
					 focal_ch
					 )
			TRG_ch = FIND_TRG.out

		} else {

			// Establish a list of TRGs using genomic phylostratigraphy (here genEra)

			// Taxids for genEra
			taxids_ch = Channel.fromPath(params.taxids)
			.splitText()
			.map{it -> it.trim()}
			.map( it -> [it.toString().split("\t")[0], it.toString().split("\t")[1]] )

			// focal taxid
			focal_taxid = taxids_ch
			.filter { name, taxid -> name == params.focal }
			.map { name, taxid -> taxid }
			.first()

			// If the user has provided a precomputed genEra output,
			if(params.genera_out){

				genera_out_ch = file(params.genera_out)

			} else {

				// A file with the paths to the neighbors CDS, plus their taxid
				neighbor_CDS_taxids = EXTRACT_CDS.out
				.map { name, fasta, gff, fai, CDS_fna, CDS_faa -> [ name, CDS_faa ] }
				.filter { name, CDS -> name != params.focal	}
				.combine(taxids_ch, by:0)
							.map { name, CDS, taxid -> "${CDS}\t${taxid}" }
							.collectFile(name: 'neighbor_CDS_taxids.txt', newLine: true)

				GENERA(
						focal_taxid,
						focal_ch
						.map{ focal_name, fasta, gff, fai, CDS_fna, CDS_faa -> CDS_faa },
						neighbor_CDS_taxids, 
						file(params.genera_db)
					  )
				genera_out_ch = GENERA.out
			}

			GENERA_FILTER(
						  genera_out_ch,
						  taxdump,
						  focal_taxid, 
						  params.trg_rank,
						  trg_node,
						  focal_mRNA_to_gene_ch
						 )
			TRG_ch = GENERA_FILTER.out
		} 


		// Then, now that we have the TRGs, find their homologs in the neighbor genomes.

		// To maximize the chances to find true homologs, all sequences (the TRG and the CDS of neighbor genomes) are elongated :
		// 100 nucleotides are added on each sides of each sequence.
		
		// The elongated parts of the TRG (of the genome of interest), are translated into the 3 possible frames to 'allow' two frameshifts in the neighbor.
		// Ad the end, 9 sequences (AA) are generated for each TRG.
		MULTIELONGATE_FOCAL_TRG( 
								TRG_ch,
								focal_ch
								)

		// Also get an elongated version of the trasnlated CDS for every genome (subjects).
		ELONGATE_CDS( EXTRACT_CDS.out )
			
		// Search for homologs of the focal's TRG among :
		// 	-the CDS (blastp),
		//	-the entiere genome (tblastn),
		// of each neighbor genome.
		// TRG collection
		BLAST	(	
				 MULTIELONGATE_FOCAL_TRG.out.TRG_multielongated_faa.first(),
				 ELONGATE_CDS.out
				)
		BLAST_out_ch = BLAST.out

	}

	// Get a list with all TRGs before additionnal filtering
	TRGS_BEFORE_STRATEGY( 
						 BLAST_out_ch.first(),
						 focal_mRNA_to_gene_ch
 						)

	// Establish the list of homologs (coding and non-coding) for each TRG and neighor genome.
	BLAST_FILTER( 
				params.focal,
				BLAST_out_ch
				)


	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Build a table to report the best hits of ech TRG in the neighbor genomes
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	// Optionnally use the checked phylogenetic tree in order to add phylogenetic distances to the table .
	if(params.tree){
		TREE_DISTANCES	(
						params.focal,
						CHECK_INPUTS.out.tree
						)
		tree_distances_ch = TREE_DISTANCES.out
	} else {
		DUMMY_DISTANCES( 
						genome_ch
						.map{ name, fasta, gff -> name }
						.collectFile(name: 'names.txt', newLine: true)
					   )
		tree_distances_ch = DUMMY_DISTANCES.out
	}

	// Build the summary table from the blast results.
	TRG_TABLE	( 
				tree_distances_ch,
				BLAST_FILTER.out.map{ name, best_hits -> best_hits }.toList()
				)
	TRG_table_ch = TRG_TABLE.out

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	If required, check the microsynteny bewteen a TRG and its non-coding homologs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	if ( params.synteny ) {

		/* 
		In order to check if a TRG and its non-coding homolog are in synteny,
		we first need to get all focal/genome orthologs pairs.
		*/

		if( params.orthodir ){
			
			//A directory with the ortholog pairs has been provided.
			// Mimic the orthologs output.
			orthologs_out_ch = EXTRACT_CDS.out
			.map { name, fasta, gff, fai, CDS_fna, CDS_faa -> [ name, file("${params.orthodir}/${params.focal}_${name}_orthologs.tsv") ] }

		} else {

			// A channel for the orthologs identification with only the name, the GFF, and the CDS_faa for each genome
			neighbors_ortho_ch = EXTRACT_CDS.out
			.map { name, fasta, gff, fai, CDS_fna, CDS_faa -> [ name, gff, CDS_faa ] }

			// A channel with only the focal genome (the GFF is not necessary here)
			focal_ortho_ch = neighbors_ortho_ch
			.map { name, gff, CDS_faa -> [ name, CDS_faa ] }
			.filter ( ~/\[${params.focal}, .*/ )
			.first() // This line converts the queue channel inton a value channel that "can be read an unlimited number of times without consuming its content" (see the doc).


			// Get all focal/genome orthologs pairs.
			// This is a subworkflow that is written in the subworkflows/orthologs.nf file.
			ORTHOLOGS(
						focal_ortho_ch,
						neighbors_ortho_ch,
					)
			orthologs_out_ch = ORTHOLOGS.out
		}

		CHECK_SYNTENY_INPUTS(
								params.focal,
								BLAST_FILTER.out
							)

		gff_ch = EXTRACT_CDS.out
		.map { name, fasta, gff, fai, CDS_fna, CDS_faa -> [ name, gff ] }

		focal_synteny_ch = gff_ch
		.filter ( ~/\[${params.focal}, .*/ )
		.first()

		synteny_main_ch = gff_ch
		.combine( orthologs_out_ch, by:0 )
		.combine( CHECK_SYNTENY_INPUTS.out, by:0 )

		// Test if synteny is conserved between the TRG genomic environment and the one of its noncoding outgroup.
		CHECK_SYNTENY(
					  focal_synteny_ch,
					  synteny_main_ch
					 )

		// Append the CHECK_SYNTENY results as a column to the TRG_table
		SYNTENY_TO_TABLE(
						 TRG_table_ch,
						 CHECK_SYNTENY.out
						 .map { out, stats -> stats }
						 .toList()
						)
		TRG_table_ch = SYNTENY_TO_TABLE.out

	}

	/*
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Interpret the results according to the given strategy and produce a final de novo genes list.
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/

	// Apply the given strategy to filter TRG_table and thus select good candidates.
	FILTER_TABLE_WITH_STRATEGY(
							   params.focal,
							   params.synteny,
							   params.strategy,
							   TRG_table_ch,
							  )

	// Remove TRGs with any isoform that has been filtered out by applying the strategy.
	FILTER_ISOFORMS(
					TRGS_BEFORE_STRATEGY.out,
					FILTER_TABLE_WITH_STRATEGY.out
				    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
