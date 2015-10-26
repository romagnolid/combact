#!/usr/bin/python 
from __future__ import print_function
import argparse
import combact
import os
import os.path
import shutil

def main():
    parser = argparse.ArgumentParser(prog="ComBact",description="Compare bacterial genomes with blast")
    parser.add_argument("-d","--database",
        nargs="+",
        help="genome(s) to be analyzed with blast")
    parser.add_argument("-o","--output",
        default="Output",
        help="output folder")
    parser.add_argument("-g","--gene-list",
        required=True,
        help="list of genes to be aligned")
    parser.add_argument("--with-silent",
        action="store_true",
        help="report silent mutations too")
    
    args = parser.parse_args()

    # output folder
    try:
        os.mkdir(args.output)
    except OSError:
        pass
    
    # make database 
    db_folder = os.path.join(args.output,"Database")
    try:    
        os.mkdir(db_folder)
    except OSError:
        pass
    
    db_file = os.path.join(db_folder,"Database.fa")
    if len(args.database) == 1: 
        shutil.copy(args.database[0], db_file)
    else:
        with open(db_file,"w") as db:
            combact.concatenate_fasta(args.database,db)
    combact.makeblastdb(db_file)
    
    # do blast search
    blast_folder = os.path.join(args.output,"Blast")
    try:    
        os.mkdir(blast_folder)
    except OSError:
        pass
    blast_file = os.path.join(blast_folder,"blast.xml")
    combact.blastn(args.gene_list, db_file, blast_file)
    
    # make reports
    reports_folder = os.path.join(args.output,"Reports")
    try:    
        os.mkdir(reports_folder)
    except OSError:
        pass
    
    full_file = open(os.path.join(reports_folder,"full.csv"),"w")
    amino_file = open(os.path.join(reports_folder,"amino.csv"),"w")
    genomes = [os.path.splitext(os.path.basename(path))[0] for path in args.database]

    if not args.report_silent:
        combact.get_mutations(blast_file, genomes, full_file, amino_file)
    else:
        combact.get_mutations(blast_file, genomes, full_file, amino_file, with_silent=True)
    
    full_file.close()
    amino_file.close()

if __name__ == "__main__":
    main()
