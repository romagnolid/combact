#!/usr/bin/env python
from __future__ import print_function
import argparse, os, shutil
from os.path import splitext, join, basename, exists
from subprocess import call
from shutil import copy
from Bio.Blast.Applications import NcbiblastnCommandline
import combact

def main():
    parser = argparse.ArgumentParser(description="ComBact 0.1")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", dest="gene_list",
        help="list of genes to be aligned")
    group.add_argument("--genbank",
        help="genbank file which genes will be extracted from")
    parser.add_argument("-d",nargs="+",dest="database",
        help="genome(s) to be analyzed with blast")
    parser.add_argument("-o", default="Output", dest="output",
        help="output folder (default=Output)")
    parser.add_argument("--with-silent", action="store_true",
        help="report silent mutations too")
    parser.add_argument("--id-cutoff", default=80, type=float,
        help="percent identity cutoff (default=80)") 
    parser.add_argument("--len-cutoff", default=70, type=float,
        help="percent length cutoff (default=70)")
    
    args = parser.parse_args()
    params = vars(args)
    
    if params.get("genbank"):
        fasta_file = splitext(basename(params["genbank"]))[0] + ".fna"
        if not exists(fasta_file):
            combact.parse_gb_file(params["genbank"], open(fasta_file,"w"))
        else:
            print(fasta_file + " already exists so it won't be overwritten. Move it, rename it or delete it.")
        gene_list = fasta_file
    else:
        gene_list = params.get("gene_list")
    
    # output folder
    try:
        os.mkdir(params["output"])
    except OSError:
        pass
    
    # make database
    db_folder = join(params["output"],"Database")
    db_file = join(db_folder,"database.fa")
    
    try:    
        os.mkdir(db_folder)
    except OSError:
        pass
    
    with open(db_file,"w") as db:
        combact.concatenate_fasta(params["database"], db)
    
    # makeblastdb
    call(["makeblastdb","-in", db_file,"-dbtype","nucl","-out", splitext(db_file)[0]])
    
    # do blast search
    blast_folder = join(params["output"],"Blast")
    blast_file = join(blast_folder,"blast.xml")
    try:    
        os.mkdir(blast_folder)
    except OSError:
        pass
    
    # blastn
    cline = NcbiblastnCommandline(query=gene_list, db=splitext(db_file)[0], outfmt=5, out=blast_file)
    print()
    print(cline)
    cline()
    
    # make reports
    reports_folder = join(params["output"],"Reports")
    try:    
        os.mkdir(reports_folder)
    except OSError:
        pass
    
    full_file = open(join(reports_folder,"full.csv"),"w")
    amino_file = open(join(reports_folder,"amino.csv"),"w")
    genomes = [splitext(basename(path))[0] for path in params["database"]]
    
    if not args.with_silent:
        combact.get_mutations(blast_file, genomes, full_file, amino_file, params["id_cutoff"], params["len_cutoff"])
    else:
        combact.get_mutations(blast_file, genomes, full_file, amino_file, params["id_cutoff"], params["len_cutoff"], with_silent=True)
    
    full_file.close()
    amino_file.close()

if __name__ == "__main__":
    main()
