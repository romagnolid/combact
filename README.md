#ComBact
ComBact (Compare Bacteria) highlights genetic differences between multiple bacterial strains of the same species (high sequences similarity is required).

## Prerequisites
[NCBI Blast+] (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
[Python] (https://www.python.org/downloads) 
[Biopython] (https://pypi.python.org/pypi/biopython)

## Installation
1. Click on the **release** tag, download latest version and unpack the tar.gz archive.

2. Add ComBact folder to your PATH with *export PATH=/path/to/ComBact:$PATH* (add the line to your **.bashrc** file to avoid typing it every time)

## Run ComBact
1. If using a reference in Genbank or Ensembl format first extract genes (and, optionally, intergenic regions) using parse\_gbk.py . Required options are (1) the reference file, and (2) the output file name; 

2. Procede to align your reference to your set of genomes using blastdb\_and\_blastn.py . As the name imply, a blast database format is created to optimize the alignment process. Required options are (1) the set of reference gene, either generated with parse\_gbk.py or provided by the user, and (2) output file name.

3. Run combact.py to produce the three tables reporting the status of each gene in every genome. Required options are (1) the alignment file produced by blastdb\_and\_blastn.py, and (2) the output folder where the tables will be created. A file listing the genomes to be analyzed is also required but being automatically produced by blastdb\_and\_blastn.py is not strictly necessary. 

## Additional utilities
1. Use combact.py with --binary flag to produce a binary matrix reporting if a gene has been found [=1] or not [=0] in a genome. Only genes present in full length (independent from SNPs and indels) are reported (fragments and multiple copies excluded).

2. Use combact.py with --fasta-core to extract sequences of genes present in every genome, representing rough versions of core genomes.

3. Length (-L)  and identity (-I) cutoffs for not reporting mutations affecting genes/fragments below those thresholds.

### Get help
A complete list with descriptions of all options can be displayed with *-h/--help* flag for all three tools.
