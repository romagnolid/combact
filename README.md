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
In order to make things easier, one can put all the input files (reference genome and genomes to be compared) into one folder and change directory accordingly.

Then, there are three possibilites:

1. To provide each input genome as a separated -i parameter (useful with few genomes):
`ComBact.sh -r ref_genome.fa -i genome1.fa -i genome2.fa`

2. To put all the input genomes into one folder and use -f parameter:
`ComBact.sh -r ref_genome.fa -f <input_folder>`

3. To provide a text file with a list of all the input genomes (one per line) and use -l parameter:
`ComBact.sh -r ref_genome.fa -l query_list.txt`

## Optional parameters
* -l parameter additionally allows users to specify a name for each input file (see query\_list.txt file for examples);

* -g parameter is an alternative to -r parameter allowing users to use a GenBank file (.gb) which is converted into a new fasta file;

* -I and -L options set two different thresholds,
percentage identity (default=80) and percentage ratio between alignment length and query length (default=70): below such values a gene is considered "Not_present";

* -s flag makes the program display silent mutations in *amino_table.csv* output file (disabled by default).

## Output
Two tab-separated files reporting gene status (WT, SNP, indel, fragment, not_present):

*full_table.csv* reports mutations at both nucleotidic and aminoacidic level (including silent mutations);

*amino_table.csv* reports only mutations at aminoacidic level (by default, silent mutations are not reported).

>If multiple copies of a gene are present, each status is reported and separated from the others by a vertical bar ("|").

>Multiple mutation events such as SNPs or indels affecting the same gene are separated by a semicolon (";").

### Example
>GeneID KlebsiellaPneumoniae\_FI KlebsiellaPneumoniae_SI

>gene1 WT WT|Fragment\_1-144|p.3A>M\_MISSENSE$p.5L_SILENT

*gene1* is wildtype in KlebsiellaPneumoniae_FI and present in three copies in KlebsiellaPneumoniae_SI: 
first copy is wildtype, second copy is a fragment and third copy is a mutant affected by two mutations (missense and silent, respectively).
