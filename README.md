# Updates
* ComBact 0.7.3 16/10/2017:
    README Updates

* ComBact 0.7.2 27/09/2016: 
    * corrected mutation report for partial alignments

---
# ComBact
ComBact (Compare Bacteria) highlights genetic differences between multiple bacterial strains of the same species (high sequences similarity is required).


ComBact is based on the idea of creating a local database using a set of genomes of interest, 
against which a set of genes is aligned 
exploiting the power and speed of the BLAST algorithm.
The references genomes must be in FASTA format (otherwise see point (1) of [Running ComBact].

## Prerequisites
* [NCBI Blast+ >= 2.5.0] (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Python] (https://www.python.org/downloads)
* [Biopython] (https://pypi.python.org/pypi/biopython)

## Installation
1. Click on the **release** tag, download latest version and unpack the *.tar.gz* archive.

2. Add ComBact folder to your PATH with `export PATH=/path/to/combact:$PATH` (add the line to your *.bashrc* file to avoid typing it every time)

## Usage
1. If using a reference in **Genbank** or **Ensembl** format, first you must extract genes 
(and, optionally, intergenic regions) using `parse_gbk.py`. 
Required options are (a) the reference file, and (b) the output file name; 

2. Align your reference to your set of genomes using `blastdb_and_blastn.py`. 
A "blast database" is created to optimize the alignment process. 
Required options are (a) the set of reference gene, 
either generated with `parse_gbk.py` or provided by the user, and (b) output file name.

3. Run `combact.py` to produce the three tables reporting the status of each gene in every genome. 
Required options are (a) the alignment file produced by `blastdb_and_blastn.py`, and (b) the output folder where the tables will be created. 
A file listing the genomes to be analyzed is also required but being automatically produced by `blastdb_and_blastn.py` is not strictly necessary.

## Options
1. Users can set length (`-L`) and identity (`-I`) thresholds to skip reporting
genes/fragments with alignment scores below given thresholds.

## Get help
A complete list with descriptions of all options can be displayed with *`-h/--help`* flag for all three tools.
