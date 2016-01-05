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
    parser.add_argument("genomes",metavar="INPUT_GENOMES",nargs="*",
        help="files that will be used to make the database")
    parser.add_argument("-l","--list",
        help="list of absolute paths of files that will be used to make the database")
    parser.add_argument("--blast-opt",type=str,default="",
        help="additional options for blastn (put them between quotes \"\")")
    args = parser.parse_args(argv)
    
    if args.genomes:
        filespaths = args.genomes
        with open("input_list.txt", "w") as inlist:
            genomes = [os.path.splitext(os.path.basename(path))[0] for path in filespaths]
            writer = csv.writer(inlist,delimiter="\t")
            for i in range(len(genomes)):
                writer.writerow([os.path.abspath(filespaths[i]), genomes[i]])
    else:
        with open(args.list) as inlist:
            reader = csv.reader(inlist,delimiter="\t")
            filespaths = []
            genomes = []
            for row in reader:
                filespaths.append(row[0])
                genomes.append(row[1])

    # create database directory
    database = os.path.join(os.getcwd(),"Database","database.fa")
    try:
        os.mkdir("Database")
    except OSError:
        pass

    # concatenate all fastas
    sequences = []
    for path in filespaths:
        genomename = os.path.splitext(os.path.basename(path))[0]
        records = SeqIO.parse(path,"fasta")
        for record in records:
            sequences.append(SeqRecord(record.seq,id=record.id,
                description=genomename))
    SeqIO.write(sequences, database,"fasta")

    # make blast database
    cline = ["makeblastdb","-dbtype","nucl","-in",database,"-out",os.path.splitext(database)[0]]
    print(" ".join(cline))
    out = subprocess.call(cline)
    if out != 0:
        sys.exit()

    # blast multifasta
    start = time.time()
    print("\nBlast alignment, current time:",time.strftime("%d/%m/%y %H:%M:%S"),file=sys.stderr)
    cline = ["blastn","-query",args.fasta,"-db",os.path.splitext(database)[0],"-outfmt","5","-out",args.output]
    print(" ".join(cline + args.blast_opt.split()))
    out = subprocess.call(cline + args.blast_opt.split())
    if out != 0:
        sys.exit()
    print("Completed in",round(time.time()-start,4),"seconds.",file=sys.stderr)

    ## blast multifasta
    #cline = NcbiblastnCommandline(query=args.fasta, db=os.path.splitext(database)[0], outfmt=5, out=args.output,num_threads=args.num_threads)
    #print("\nBlast alignment, current time:",time.strftime("%d/%m/%y %H:%M:%S"),file=sys.stderr)
    #start = time.time()
    #print(str(cline))
    #cline()
    #print("Completed in",round(time.time()-start,4),"seconds.",file=sys.stderr)
            
if __name__ == "__main__":
    main(sys.argv[1:])
