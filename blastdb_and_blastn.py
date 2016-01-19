#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import os
import os.path
import subprocess
import time
import sys
import csv

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta",metavar="INPUT_FILE",
        help="multifasta with genes of interest")
    parser.add_argument("output",metavar="OUTPUT_FILE",
        help="blast output filename (format xml)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g","--genomes",metavar="INPUT_GENOMES",nargs="+",
        help="files that will be used to make the database")
    group.add_argument("-i","--inlist",metavar="INPUT_LIST",
        help="list of absolute paths of files that will be used to make the database")
    parser.add_argument("--blastn",type=str,default="",
        help="additional options for blastn (between quotes \"\")")
    parser.add_argument("--makeblastdb",type=str,default="",
        help="additional options for blastn (between quotes \"\")")
    args = parser.parse_args(argv)

    if args.genomes:
        # create inlist to be used later
        paths = args.genomes
        names = [os.path.splitext(os.path.basename(path))[0] for path in paths]
        genomes = zip(paths, names)
        with open("inlist.txt", "w") as infile:
            infile.writelines(["{str}\n".format(str=name) for name in names])
    elif args.inlist:
        with open(args.inlist) as infile:
            paths = [line.split()[0] for line in infile]
            names = [os.path.splitext(os.path.basename(path))[0] for path in paths]
            genomes = zip(paths,names)

    # create database directory
    database = os.path.join(os.getcwd(),"Database","database.fa")
    try:
        os.mkdir("Database")
    except OSError:
        pass

    # concatenate all fastas
    sequences = []
    for path in paths:
        genomename = os.path.splitext(os.path.basename(path))[0]
        records = SeqIO.parse(path,"fasta")
        for record in records:
            sequences.append(SeqRecord(record.seq,id=record.id,
                description=genomename))
    SeqIO.write(sequences, database,"fasta")

    # make blast database
    cline = ["makeblastdb","-dbtype","nucl","-in",database,"-out",os.path.splitext(database)[0]]
    print(" ".join(cline + args.makeblastdb.split()))
    out = subprocess.call(cline + args.makeblastdb.split())
    if out != 0:
        sys.exit()

    # blast multifasta
    start = time.time()
    print("\nBlast alignment, current time:",time.strftime("%d/%m/%y %H:%M:%S"),file=sys.stderr)
    cline = ["blastn","-query",args.fasta,"-db",os.path.splitext(database)[0],"-outfmt","5","-out",args.output]
    print(" ".join(cline + args.blastn.split()))
    out = subprocess.call(cline + args.blastn.split())
    if out != 0:
        sys.exit()
    print("Completed in",round(time.time()-start,4),"seconds.",file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])
