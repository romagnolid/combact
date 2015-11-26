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
    parser.add_argument("--add-igr", action="store_true",
        help="extract intergenic region from genbank file")
    parser.add_argument("--add-silent", action="store_true",
        help="report silent mutations too")
    parser.add_argument("--id-cutoff", default=80, type=float,
        help="percent identity cutoff (default=80)") 
    parser.add_argument("--len-cutoff", default=70, type=float,
        help="percent length cutoff (default=70)")
    
    args = parser.parse_args()
    params = vars(args)
    
    # output folder
    try:
        os.mkdir(params["output"])
    except OSError:
        pass
    
    # genbank or gene list
    if params.get("genbank"):
        fasta_file = join(params["output"],splitext(basename(params["genbank"]))[0] + ".fna")
        if not exists(fasta_file):
            if params["add_igr"]:
                combact.parse_gb_file(params["genbank"], open(fasta_file,"w"), True)
            else:
                combact.parse_gb_file(params["genbank"], open(fasta_file,"w"))
        gene_list = fasta_file
    else:
        gene_list = params.get("gene_list")
    
    # make database
    db_folder = join(params["output"],"Database")
    db_file = join(db_folder,"database.fa")
    
    try:    
        os.mkdir(db_folder)
    except OSError:
        pass
    
    if not exists(db_file):
        combact.concatenate_fasta(params["database"], open(db_file,"w"))
        # makeblastdb
        call(["makeblastdb","-in", db_file,"-dbtype","nucl","-out", splitext(db_file)[0]])
    else:
        print("\nA database already exists. Delete the old file to make a new one.")
    
    # do blast search
    blast_folder = join(params["output"],"Blast")
    blast_file = join(blast_folder,"blast.xml")
    try:    
        os.mkdir(blast_folder)
    except OSError:
        pass
    
    # blastn
    if not exists(blast_file):
        cline = NcbiblastnCommandline(query=gene_list, db=splitext(db_file)[0], outfmt=5, out=blast_file)
        print("\n" + cline)
        cline()
    else:
        print("\nBlast output already exists. Delete the old file to make a new one.")
    
    # make reports
    reports_folder = join(params["output"],"Reports")
    try:    
        os.mkdir(reports_folder)
    except OSError:
        pass
    
    full_file = open(join(reports_folder,"full.csv"),"w")
    amino_file = open(join(reports_folder,"amino.csv"),"w")
    genomes = [splitext(basename(path))[0] for path in params["database"]]
    
    if not args.add_silent:
        combact.get_mutations(blast_file, genomes, full_file, amino_file, params["id_cutoff"], params["len_cutoff"])
    else:
        combact.get_mutations(blast_file, genomes, full_file, amino_file, params["id_cutoff"], params["len_cutoff"], add_silent=True)
    
    full_file.close()
    amino_file.close()

if __name__ == "__main__":
    main()
