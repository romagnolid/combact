## ComBact 0.8
ComBact (Compare Bacteria) is a set of python scripts to highlight genetic
differences between multiple bacterial strains of the same species (high
sequences similarity is required).

ComBact is based on the idea of creating a local database using a set of
genomes of interest against which a set of genes is aligned, exploiting the
accuracy and speed of the BLAST algorithm.  The reference genome must be in
FASTA format (otherwise see point (1) of [Usage]).

####
Mutations nomenclature is loosely based on HGVS raccomendations
(http://varnomen.hgvs.org/)

## Prerequisites
* [NCBI Blast+ >= 2.7.0](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Python2](https://www.python.org/downloads)
* [Biopython](https://pypi.python.org/pypi/biopython)

## Installation
1. Clone the package or download the latest release and unpack the *.tar.gz*
archive.

2. Add ComBact folder to your PATH with `export PATH=/path/to/combact:$PATH`
(add the line to your *.bashrc* file to avoid typing it every time).

## Usage
1. If using a reference in **Genbank** or **Ensembl** format, first you must
extract genes (and, optionally, intergenic regions) using `parse_gbk.py`.
Required options are the reference file, and the output file name;

2. Align your reference to your set of genomes using `blastdb_and_blastn.py`.
First a "blast database" is created to optimize the alignment process. Then the
alignment process starts and the results are saved in two files: an ASN.1 blast
format, from which users can create their preferred output format using
`blast_formatter`, and an XML output which is then parsed by the third script.
Required options are the fasta file with the reference genes, the name of the
output file (without any extension) and the text file reporting the absolute
paths to each genome to analyze. The latter will be used in the next step too.

3. Run `combact.py` to produce the three tables reporting the mutation status
of each gene in the selected genomes.  Required options are the XML file
produced by `blastdb_and_blastn.py`, the output folder where the output tables
are created, and the text file listing the genomes to be analyzed.
Users can provide subsets of such list to generate tables with only a selection
of genomes.

## Options
Users can set length (`-L`) and identity (`-I`) thresholds to
avoid reporting genes/fragments with alignment scores below given thresholds.

## Example of mutations
- SNP: c.[64C>T;102T>G]
- INDEL: c.[64_66delG;102T>G]
- FRAGMENT: c.[64C>T;102T>G](25:104)
