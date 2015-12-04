#ComBact README
ComBact (Compare Bacterial genomes) highlights genetic differences between multiple bacterial strains of the same species (high sequences similarity is required).

## Requirements
The following software are required:

* [Python] (https://www.python.org/downloads)
* [Biopython] (https://pypi.python.org/pypi/biopython)
* [NCBI Blast+] (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## Installation
1. Download and decompress the zip archive.

2. Include folder path to your $PATH variable using
`export PATH=/path/to/ComBact:$PATH`

## Usage
1. List every files that will make the database with the -d option (optionally using Unix-like glob-patterns):
`run_combact.py -d ecoli1.fa ecoli2.fa ecoli3.fa ...`
`run_combact.py -d ECOLI_folder/* ...`

2. Provide a text file with a list of the genomes (one per line) with -l option:
`run_combact.py -l query_list.txt ...`

## Other parameters
* -r | --reference: reference genome (fasta format)

* -g | --genbank: reference genome (Genbank format: will be converted into fasta format);

* -I and -L options: percentage identity cutoff (default=80) and alignment and query length percentage ratio (default=70) (genes below cutoff will be considered "absent");

* --add-silent: report silent mutations (beta-version; affect only "full" and "amino" reports);

* --add-igr: either extract intergenic regions from genbank file or not.

## Output
Three tab-separated files reporting gene status (wt, SNP, indel, wt\_fragment, fragment\_mut).

Two text files reporting genes either being wild-type in all genomes or absent (according to BLAST) in all genomes.

*full.csv* reports mutations at both nucleic and aminoacid level;

*nucl.csv* reports mutations at nucleic level;

*amino.csv* reports mutations at aminoacid level.

If multiple copies of a gene are present, each status is reported and separated from the others by "|".

Multiple mutation events such as SNPs or indels affecting the same gene are separated by a semicolon (";").
