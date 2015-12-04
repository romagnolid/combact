#!/usr/bin/env python
from __future__ import print_function
import argparse, os, shutil
import time, datetime
from os.path import splitext, join, basename, exists
from subprocess import call
from shutil import copy
from Bio.Blast.Applications import NcbiblastnCommandline
import combact

def main(argv=None):
    parser = argparse.ArgumentParser(description="ComBact_0.3")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-r", "--reference", dest="gene_list",
        help="reference genome")
    group.add_argument("-g", "--genbank",
        help="reference genome (Genbank format)")

    parser.add_argument("-l",nargs=1,dest="file_list",
        help="genome(s) to be analyzed with blast")

    parser.add_argument("-d",nargs="+",dest="database",
        help="genome(s) to be analyzed with blast")
    parser.add_argument("-o", default="Combact_output", dest="output",
        help="output folder (default=Combact_output)")
    parser.add_argument("-L","--length-cutoff", default=70, type=float,
        help="percent length cutoff (default=70)")
    parser.add_argument("-I","--identity-cutoff", default=80, type=float,
        help="percent identity cutoff (default=80)") 
    parser.add_argument("--add-igr", action="store_true",
        help="extract intergenic region from genbank file")
    parser.add_argument("--add-silent", action="store_true", dest="silent_opt",
        help="report silent mutations too")
    
    if argv is not None:
        args = parser.parse_args(argv.split())
    else:
        args = parser.parse_args()

    params = vars(args)

#    print(params)
#def nope():

    print("\nAnalysis started at " + time.ctime())
    start = time.time()

    # output folder
    try:
        os.mkdir(params["output"])
    except OSError:
        pass
    
    # genbank or gene list
    if params.get("genbank"):
        fasta_file = join(params["output"], splitext(basename(params["genbank"]))[0] + ".fna")
        if not exists(fasta_file):
            print("Parsing Genbank file...")
            if params["add_igr"]:
                combact.parse_gb_file(params["genbank"], open(fasta_file,"w"), True)
            else:
                combact.parse_gb_file(params["genbank"], open(fasta_file,"w"))
            print("Done", str(datetime.timedelta(seconds=time.time()-start)))
        else:
            print("\nWarning: Genbank file already been parsed. Delete it to create a new one.")
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
    
    if params["database"]:
        files = params["database"]
    elif params["file_list"]:
        with open(params["file_list"][0]) as fl:
            lines = fl.readlines()
            files = [line.split()[0] for line in lines]
    # makeblastdb
    if not exists(db_file):
        combact.cat_fasta(files, open(db_file,"w"))
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
        print("Done", str(datetime.timedelta(seconds=time.time()-start)))

        print("\nCreating xml blast output...")
        call(["blast_formatter", "-archive", blastAsn_file, "-out", blastXml_file, "-outfmt", "5"])
        print("Done", str(datetime.timedelta(seconds=time.time()-start)))

        print("\nCreating text blast output...")
        call(["blast_formatter", "-archive", blastAsn_file, "-out", blastTxt_file, "-outfmt", "0"])
        print("Done", str(datetime.timedelta(seconds=time.time()-start)))

    else:
        print("\nWarning: blast output already exist. Delete it to create a new one.")
    
    # make reports
    reports_folder = join(params["output"],"Reports")
    try:    
        os.mkdir(reports_folder)
    except OSError:
        pass
    
    print("\nParsing blast.xml output...")
    genomes = [splitext(basename(path))[0] for path in files]
    combact.get_mutations(blastXml_file, genomes, reports_folder, params["identity_cutoff"], params["length_cutoff"], params["silent_opt"])
    
    print("Done", str(datetime.timedelta(seconds=time.time()-start)))

if __name__ == "__main__":
    main()
