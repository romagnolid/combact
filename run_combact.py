#!/usr/bin/env python
from __future__ import print_function
import argparse, os, shutil
from os.path import splitext, join, basename, exists
from subprocess import call
from shutil import copy
from Bio.Blast.Applications import NcbiblastnCommandline
import combact

def main():
    parser = argparse.ArgumentParser(description="ComBact 0.2")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", dest="gene_list",
        help="list of genes to be aligned")
    group.add_argument("-g", "--genbank",
        help="genbank file which genes will be extracted from")
    parser.add_argument("-d",nargs="+",dest="database",
        help="genome(s) to be analyzed with blast")
    parser.add_argument("-o", default="Combact_output", dest="output",
        help="output folder (default=Combact_output)")
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
        print("Parsing Genbank file...")
        if not exists(fasta_file):
            if params["add_igr"]:
                combact.parse_gb_file(params["genbank"], open(fasta_file,"w"), True)
            else:
                combact.parse_gb_file(params["genbank"], open(fasta_file,"w"))
        gene_list = fasta_file
        print("Done.")
    else:
        gene_list = params.get("gene_list")
    
    # make database
    db_folder = join(params["output"],"Database")
    db_file = join(db_folder,"database.fa")
    try:    
        os.mkdir(db_folder)
    except OSError:
        pass
    
    # makeblastdb
    if not exists(db_file):
        combact.concatenate_fasta(params["database"], open(db_file,"w"))
        call(["makeblastdb","-in", db_file,"-dbtype","nucl","-out", splitext(db_file)[0]])
    else:
        print("\nWarning: database already exists. Delete it to create a new one.")
    
    # do blast search
    blast_folder = join(params["output"],"Blast")
    blastAsn_file = join(blast_folder,"blast.asn")
    blastXml_file = join(blast_folder,"blast.xml")
    blastTxt_file = join(blast_folder,"blast.txt")
    try:    
        os.mkdir(blast_folder)
    except OSError:
        pass
    
    if not exists(blastAsn_file):
        # blastn
        cline = NcbiblastnCommandline(query=gene_list, db=splitext(db_file)[0], outfmt=11, out=blastAsn_file)
        print("\nBlast command line:\n" + str(cline))
        cline()
        print("\nCreating xml blast output.")
        call(["blast_formatter", "-archive", blastAsn_file, "-out", blastXml_file, "-outfmt", "5"])
        print("\nCreating text blast output.")
        call(["blast_formatter", "-archive", blastAsn_file, "-out", blastTxt_file, "-outfmt", "0"])
    else:
        print("\nWarning: blast output already exist. Delete it to create a new one.")
    
    # make reports
    reports_folder = join(params["output"],"Reports")
    try:    
        os.mkdir(reports_folder)
    except OSError:
        pass
    
    full_file = open(join(reports_folder,"full.csv"),"w")
    nucl_file = open(join(reports_folder,"nucl.csv"),"w")
    amino_file = open(join(reports_folder,"amino.csv"),"w")
    genomes = [splitext(basename(path))[0] for path in params["database"]]
    
    print("\nParsing blast.xml output...")
    if not params["add_silent"]:
        combact.get_mutations(blastXml_file, genomes, full_file, nucl_file, amino_file, params["id_cutoff"], params["len_cutoff"])
    else:
        combact.get_mutations(blastXml_file, genomes, full_file, nucl_file, amino_file, params["id_cutoff"], params["len_cutoff"], add_silent=True)
    
    full_file.close()
    nucl_file.close()
    amino_file.close()
    print("Done.")

if __name__ == "__main__":
    main()
