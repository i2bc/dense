## Input/output options
Define where the pipeline should find input data and save output data.

| Parameter | Description | Default |
| --- | --- | --- |
| outdir | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. |  |
| gendir | The input directory that contains a genomic FASTA file ('.fna','.fasta') and a GFF3 annotation file ('.gff','.gff3') for each genome (focal and neighbors). |  |

## Required input values


| Parameter | Description | Default |
| --- | --- | --- |
| focal | The name of the focal genome (the one whose CDS will be tested). |  |
| strategy | The strategy number to apply to detect denovo genes :  1-> has a non-coding hit 'outgroup' (see https://github.com/Proginski/dense/tree/master#outgroup) i.e. has a non-coding hit in a genome without any CDS hit AND that genome has an older MRCA than any genome with a CDS hit ;  2-> has a non-coding hit, i.e. TRG has a non-coding hit in a genome without any CDS hit ;  3-> has no (but self) coding hit and has a non-coding hit | 1 |

## Other important parameters


| Parameter | Description | Default |
| --- | --- | --- |
| tree | The phylogenetic tree that shows relations between the genomes (Newick format). E.g. "((Ptep:0.75680000,Pruf:0.75680000)'210':0.38205000,Pfoa:1.13885000)'220';" |  |
| num_outgroups | The required number of outgroups with a non-coding match for the strategies 1 and 3 (see https://github.com/Proginski/dense/tree/master#outgroup) | 1 |
| help | Throws this page. |  |

## Institutional config options
Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Default |
| --- | --- | --- |
| config_profile_name | Institutional config name. |  |
| config_profile_description | Institutional config description. |  |

## Max job request options
Set the top limit for requested resources for any single job.

| Parameter | Description | Default |
| --- | --- | --- |
| max_cpus | Maximum number of CPUs that can be requested for any single job. | 16 |
| max_memory | Maximum amount of memory that can be requested for any single job. | 128.GB |
| max_time | Maximum amount of time that can be requested for any single job. | 240.h |

## TRG list options


| Parameter | Description | Default |
| --- | --- | --- |
| trgsblastdir | A directory with the (precomputed) BLAST output necessary for TRG homologs detection. Two files per neighbor genome, must be named 'TRG_multielongated_blastp_${neighbor}_CDS_elongated.out' and 'TRG_multielongated_tblastn_${neighbor}_genome.out'. Incompatible with the following parameters. |  |
| trgs | A text file with a predefined list of CDS to consider as TRGs (incompatible with the following parameters). |  |
| genera_out | A '.tsv' file with precomputed gene ages from genEra. Makes '--genera_db' useless. |  |
| genera_db | Path to the diamond db with taxonomy (e.g. '../DBs/nr.dmnd'). Necessary if a list of TRGs is not provided ('--TRGs'). |  |
| genera_fast | Use a custom version of GenEra to deal with large genomes. Warning : may open as many files as the number of CDSs. | False |
| taxdump | The taxdump directory path (otherwise downloaded). |  |
| taxids | A '.tsv' file with two columns : col1 = genome name, col2 = taxid. Must include all genomes (focal and neighbors). |  |
| trg_node | A taxonomic node (e.g. Mammalia) to filter CDS into TRGs. CDS associated with this node or on of ots children will be considered as TRGs (incompatible with '--trg_rank'). |  |
| trg_rank | A taxonomic rank (e.g. 'order') to filter CDS into TRGs. CDS associated with this node or on of ots children will be considered as TRGs (incompatible with '--trg_node'). | genus |

## Synteny checking


| Parameter | Description | Default |
| --- | --- | --- |
| synteny | Whether or not the check is TRGs are in synteny with their non-coding homolog(s). Required by the following parameters. | True |
| synteny_window | The number of flanking genes to collect on each side (5' and 3') of the TRG and on each side of its non-coding hit. | 4 |
| synteny_anchors | The minimum number of genes on each side of the TRG that should have their ortholog on one side of its non-coding hit. E.g. : out of 4 flanking genes (set with synteny_window), 3 upstream genes have their ortholog among the flanking genes of the non-coding hit, and 2 for downtream genes. The synteny will be assumed if synteny_anchors is set to 1 or 2, but not 3 or more, because there are only 2 'anchors' on the downstream side of the query. | 1 |
| orthodir | A directory with the precomputed pairs of orthologous genes for each focal-neighbor genomes pair. One file per genomes pair, must be named '${focal}_${neighbor}_orthologs.tsv'. Incompatible with the following parameters. |  |
| orthotool | The tool to perform alignments for orthologs definition. | diamond |
| diamond_sens | Diamond sensitivity. |  |

