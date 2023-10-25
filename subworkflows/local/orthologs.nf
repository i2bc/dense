/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BLAST                } from '../../modules/local/orthologs_modules.nf'
include { DIAMOND_BLAST        } from '../../modules/local/orthologs_modules.nf'
include { BEST_HITS            } from '../../modules/local/orthologs_modules.nf'
include { MRNA_TO_GENE         } from '../../modules/local/orthologs_modules.nf'
include { RECIPROCAL_BEST_HITS } from '../../modules/local/orthologs_modules.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow ORTHOLOGS {

	/* 
	This workflow finds the orthologs for a genome and its neighbors.
	It takes two tuple as input (see 'take') and outputs a tuple per neighbor on the model : 
	[ neighbor_name, pairs_of_orthologs.tsv ].
	*/

	take :
		focal_ortho_ch 		// tuple [name, CDS.faa ]
		neighbors_ortho_ch	// tuple [name, gff, CDS.faa ]
	
	main :


	// First, get the pairwise alignments.
	if ( params.blastdir ) {
		
		//A directory with the BLAST outputs was provided.
		// Mimic the BLAST process output.
		BLAST_out_ch = neighbors_ortho_ch
		.map { name, gff, CDS_faa -> [ params.focal, name, file("${params.blastdir}/${params.focal}_CDS_BLASTp_${name}_CDS.out"), file("${params.blastdir}/${name}_CDS_BLASTp_${params.focal}_CDS.out") ] }

	} else {

		// For each focal genome/genome pair, perform two BLASTp :
		// - focal genome translated CDS against genome translated CDS
		// - genome translated CDS against focal genome translated CDS
		neighbor_blast_ch = neighbors_ortho_ch
		.map { name, gff, CDS_faa -> [ name, CDS_faa ] }

		if ( params.blasttool == "diamond" ){

			if (params.diamond_sens) {
				sensitivity = "--${params.diamond_sens}"
			} else {
				sensitivity = ""
			}
			DIAMOND_BLAST(
							focal_ortho_ch,
							neighbor_blast_ch,
							sensitivity
						  )
			BLAST_out_ch = DIAMOND_BLAST.out

		}
		if ( params.blasttool == "blast" ){

			BLAST(
					focal_ortho_ch,
					neighbor_blast_ch
				 )
			BLAST_out_ch = BLAST.out

		}

	}
	
	
	// For each blast output file get a two columns tsv file with each query and its best hit's subject.
	BEST_HITS( BLAST_out_ch )


	// For each GFF3 file, get a CDS-to-gene mapping (two cols tsv file). 
	// It will be used to "pool" different mRNA from the same locus (gene) together in the orthologs definition.
	gff_ch = neighbors_ortho_ch
	.map { name, gff, CDS_faa -> [ name, gff ] }

	MRNA_TO_GENE( gff_ch )
	

	// For each focal-genome/neighbor pair, get a list of the reciprocal gene best hits (i.e. the list of orthologs).
	focal_mRNA_to_gene_ch = MRNA_TO_GENE.out
	.filter ( ~/\[${params.focal}, .*/ )
	.first()

	final_input_ch = BEST_HITS.out.mainout
  	// Combine according to a key that is the first value of every first element (here 'name').
  	.combine( MRNA_TO_GENE.out, by: 0 )


	/*
	For a genome F and its neighbor N we have two .tsv files :
	- first .tsv  : F_CDS (col 1) and N_CDS_best_hit (col2)
	- second .tsv : N_CDS (col 1) and F_CDS_best_hit (col2)
	Turn these into (virtually):
	- first .tsv  : F_gene (col 1) and N_gene_best_hit (col 2)
	- second .tsv : N_gene (col 1) and F_gene_best_hit (col 2)
	Keep the common pairs between the two files (.i.e : the reciprocal best hits).
	*/
	RECIPROCAL_BEST_HITS( 
							focal_mRNA_to_gene_ch,
							final_input_ch
						)

	emit :
		RECIPROCAL_BEST_HITS.out
	
}
