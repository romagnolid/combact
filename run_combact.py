#!/usr/bin/env python
from __future__ import print_function
import argparse
import combact
import os
import os.path
import shutil

def main():
    parser = argparse.ArgumentParser(description="Compare bacterial genomes with blast")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gene-list",
        help="list of genes to be aligned")
    group.add_argument("--gb",
        help="genbank file which genes will be extracted from")
    parser.add_argument("-d","--database",
        nargs="+",
        help="genome(s) to be analyzed with blast")
    parser.add_argument("-o","--output",
        default="Output",
        help="output folder")
    parser.add_argument("--with-silent",
        action="store_true",
        help="report silent mutations too")
    
    args = parser.parse_args()

    if args.gb:
        fasta_file = os.path.splitext(args.gb)[0] + ".fna"
        if os.path.exists(fasta_file):
            raise IOError("Trying to create " + fasta_file + "but already present.")
        else:
            combact.parse_gb_file(args.gb, open(fasta_file,"w"))
        gene_list = fasta_file

    else:
        gene_list = args.gene_list

    # output folder
    try:
        os.mkdir(args.output)
    except OSError:
        pass
    
    # make database 
    db_folder = os.path.join(args.output,"Database")
    db_file = os.path.join(db_folder,"Database.fa")
    try:    
        os.mkdir(db_folder)
        combact.concatenate_fasta(args.database,open(db_file,"w"))
        combact.makeblastdb(db_file)

    except OSError:
        answer = raw_input("Database already exists. Overwrite it? [y/n]: ")
        if answer and answer.upper().startswith("Y"):
            if len(args.database) == 1:
                shutil.copy(args.database[0], db_file)
            else:
                combact.concatenate_fasta(args.database,open(db_file,"w"))
                combact.makeblastdb(db_file)
    
    # do blast search
    blast_folder = os.path.join(args.output,"Blast")
    try:    
        os.mkdir(blast_folder)
    except OSError:
        pass
    blast_file = os.path.join(blast_folder,"blast.xml")
    combact.blastn(gene_list, db_file, blast_file)
    
    # make reports
    reports_folder = os.path.join(args.output,"Reports")
    try:    
        os.mkdir(reports_folder)
    except OSError:
        pass
    
    full_file = open(os.path.join(reports_folder,"full.csv"),"w")
    amino_file = open(os.path.join(reports_folder,"amino.csv"),"w")
    genomes = [os.path.splitext(os.path.basename(path))[0] for path in args.database]

    if not args.with_silent:
        combact.get_mutations(blast_file, genomes, full_file, amino_file)
    else:
        combact.get_mutations(blast_file, genomes, full_file, amino_file, with_silent=True)
    
    full_file.close()
    amino_file.close()

if __name__ == "__main__":
    main()
