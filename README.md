# ![DENSE](docs/images/Dense_logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**WARNING : THIS REPOSITORY IS under CONSTRUCTION**
** BOTH DOCUMENTATION AND SCRIPTS ARE NOT COMPLETED OR UP TO DATE**
** PLEASE WAIT FOR A FIRST RELEASE | EARLY 2024 **

## Introduction

**DENSE** is a pipeline that detects genes that have emerged *de novo* (from non-coding DNA regions), based on phylostratigraphy and synteny.

![dag.png](docs/images/flowchart_vs6.png)

**DENSE** uses a genome of interest (focal) and its phylogenetic neighbors (genomes FASTA and GFF3 annotation files). The pipeline includes the following steps :

* **A :** extracts the coding sequences of all protein coding genes in the focal genome and search for homologs among the Refseq Non-redundant protein database (NR), and the neighbor genomes.
  
* **B :** based on the previous step, selects genes that are taxonomically restricted (TRG)
  
* **C :** selects TRG with homology in non-coding regions of neighbor genomes. If a phylogenetic tree was provided, DENSE can require from these genomes to be 'outgroup', meaning that there are more distant from the focal genome that any neighbor actually sharing the gene.
  
* **D :** DENSE finally determines whether the homologous non-coding regions are in synteny with their TRG (the step can be switch off).  
It generates a file containing all the genes that have emerged *de novo*.

## Table of contents

<!--ts-->
- [](#)
  - [Introduction](#introduction)
  - [Table of contents](#table-of-contents)
  - [Set-up](#set-up)
    - [1. Nextflow](#1-nextflow)
    - [2. Container manager](#2-container-manager)
    - [3. Download the NR (not mandatory)](#3-download-the-nr-not-mandatory)
  - [Input files](#input-files)
  - [Usage](#usage)
    - [Lucy example](#lucy-example)
      - [command](#command)
      - [config file](#config-file)
    - [Luca example](#luca-example)
      - [command](#command-1)
      - [config file](#config-file-1)
  - [Options](#options)
  - [Pipeline output](#pipeline-output)
  - [Credits](#credits)
  - [Citations](#citations)
<!--te-->

## Set-up

### 1. Nextflow

Before anything, you need to have an recent Nextflow installed.
> If you do not have Nextflow yet, you can find simple instructions here : [this page](https://www.nextflow.io/docs/latest/getstarted.html).  

In order to use the latest Nextflow version, you should use:
```bash
nextflow self-update
```
> [!IMPORTANT] 
> DENSE **requires** Nextflow >=23.04.3. A previous version could lead to errors.

To test your Nextflow installation you can use : 
```bash
nextflow run hello
```

### 2. Container manager

In order to use DENSE in an fully-ready and reproducible environment, you need to have a container manager installed on your machine.  
You can use any of the following :
* [Apptainer](https://apptainer.org/docs/user/latest/quick_start.html)
* Singularity
* [Docker](https://www.docker.com/get-started/)

You can now test **DENSE** on the example data with the following command :
```bash
nextflow run proginski/dense -profile <DOCKER|APPTAINER|SINGULARITY>,test
```

### 3. Download the NR (not mandatory)

In order to detect taxonomically restricted genes (TRG), DENSE uses [GenEra](https://github.com/josuebarrera/GenEra) to search the Refseq Non-redundant protein database (NR). 

To download and properly install the NR along with taxonomic data, you can follow [these instructions](https://github.com/josuebarrera/GenEra/wiki/Setting-up-the-database(s)).

>The downloading step can take a couple of hours, but is necessary to assess the absence of homology of your genes candidate to any other known protein coding gene.  
>You can ignore this step if you want to use you own user-defined TRG list instead (see Usage).

## Input files
To run DENSE you always need a directory that contains a genomic FASTA file ('.fna','.fasta') and a GFF3 annotation file ('.gff','.gff3') for each genome (focal and neighbors, e.g. : the mouse and some close rodents). `--gendir`
>GFF3 files must have a classical CDS < mRNA < gene parent relationship between features.

If you want to use DENSE the most complete way, you also need :
* The phylogenetic tree that shows relations between the genomes (Newick format) `--tree`
* A '.tsv' file with two columns : col1 = genome name, col2 = taxid. Must include all genomes (focal and neighbors) `--taxids`

> Get the taxid of your species:
> * GFF3 files from NCBI have a header line.
> ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=\<TAXID\>
> * Find your organism on the [NCBI Taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).
>
>Then you can write a file similar to this one : 
>```
>Droso_melanogaster  7227
>Droso_virilis   7244
>Droso_simulans  7240
>...
>```

## Usage
### Lucy example
Lucy has a favorite species. She wants to collect genes from that species with the best indications of *de novo* emergence.  
Therefore, she runs a complete DENSE analysis.  
Since her HPC's admin does not like Docker (they *all* do), she uses [Apptainer](https://apptainer.org/docs/user/latest/quick_start.html).
#### command
```
nextflow run proginski/dense -profile apptainer -c Lucy.config
```
#### config file
Lucy.config content : 
```
params {
    
    gendir    = "../GENOMES/"     // a directory that contains favorite.fna, favorite.gff3, cousin1.fna, cousin1.gff3, cousin2.fna, cousin2.gff3, ...
    focal     = "favorite"       // the name of the focal genome

    tree      = "family_tree.nwk" // a tree with the same names as the genome files 
    genera_db = "../../../nr/"    // the path to the 'nr.dmd' parent directory
    taxids    = "taxids.tsv"      // see the Input files section
    
}
```
### Luca example
Luca has dozen of annotated strains from its most cherished Yeast. 
He wants to know if its first strain has genes that seem to have emerged de novo by comparison with the eleven other strains.  
He already has a list of orphan genes for this yeast, and so he provides it to DENSE (basically skip step A and B) (`trgs`).  
He does not know the evolutionary relationship between the genomes (no `tree` and `strategy = 2`).  
He does not care about checking the synteny (`synteny = false`).  
He changed his mind about this options in the middle of a first analysis, so this time he use `-resume` to reuse pre-computed steps.
#### command
```
nextflow run proginski/dense -profile docker -c Luca.config -resume
```
#### config file
Luca.config content :
```
params {
    
    gendir   = "input_file/"              // a directory that contains strain1.fna, strain1.gff3, strain2.fna, strain2.gff3, etc...
    focal    = "strain1"                  // the name of the focal genome

    trgs     = "list_of_orphan_genes.txt" // see the Input files section
    strategy = 2                          // select any TRG with a non-coding homolog region (and no coding homolog) of a neighbor.
    synteny  = false                      // turn off synteny checking 
    
}

```
> Find out more ways to use options in Nextflow : [configs](https://www.nextflow.io/docs/latest/config.html) 

## Options
see [Options](nextflow_schema.json)

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/DENSE/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/DENSE/output).

## Credits

DENSE was originally written by Paul Roginski.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Citations

If you use  DENSE for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.