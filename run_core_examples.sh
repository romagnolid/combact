## start extracting selected loci from genbank
parse_gbk.py PAO1.gb -o PAO1_core.fa -t tag_list.txt

## create blasta database with genomes of interest and align core genes
makeblastdb_and_blastn.py core_genes.fa subset/*

## parse blast xml results to produce core_genomes
core_genomes.py blast.xml CoreGenomes
