## start extracting loci from genbank

## create blasta database with genomes of interest and align genes of interest
makeblastdb_and_blastn.py core_genes.fa examples/subset

## parse blast xml results to compare genomes
combact.py blast.xml
