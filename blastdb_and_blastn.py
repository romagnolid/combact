#!/usr/bin/python

from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import os.path
import subprocess
import time
import sys
import csv

def main(argv=None):
    parser = argparse.ArgumentParser(description="ComBact-0.8")
    parser.add_argument("fasta", metavar="INPUT_FILE",
        help="multifasta with genes of interest")
    parser.add_argument("output", metavar="OUTPUT_FILE",
        help="blast output filename (with out extention)")
    parser.add_argument("-i", "--input_list", metavar="INPUT_LIST", required = True, 
        help="list of absolute file paths of genomes to be used as database")
    parser.add_argument("-d", "--database", metavar="DATABASE_DIRECTORY",
        default = "database",
        help="directory where database will be stored (default = 'database'")
    parser.add_argument("--makeblastdb", type=str, default="",
        help="additional options for makeblastdb (between quotes \"\")")
    parser.add_argument("--blastn", type=str, default="",
        help="additional options for blastn (between quotes \"\")")
    args = parser.parse_args(argv)

    with open(args.input_list) as infile:
        paths = [line.split()[0] for line in infile]
        names = [os.path.splitext(os.path.basename(path))[0] for path in paths]
        genomes = zip(paths, names)

    # create database directory
    database = os.path.join(os.getcwd(), args.database, 
        os.path.basename(args.database) + ".fa")
    if not os.path.isdir(args.database):
        os.mkdir(args.database)

    # concatenate all fastas
    sequences = []
    for path in paths:
        records = SeqIO.parse(path, "fasta")
        for record in records:
            sequences.append(SeqRecord(record.seq, id=record.description,
            description=""))
    SeqIO.write(sequences, database, "fasta")

    # make blast database
    cline = ["makeblastdb", "-parse_seqids", "-dbtype", "nucl", "-in", database,
        "-out", os.path.splitext(database)[0]]
    print(" ".join(cline + args.makeblastdb.split()))
    out = subprocess.call(cline + args.makeblastdb.split())

    # blast multifasta
    start = time.time()
    print("\nBlast alignment, current time:", time.strftime("%d/%m/%y %H:%M:%S"))

    # output format asn
    cline = ["blastn", "-query", args.fasta, "-db", os.path.splitext(database)[0],
        "-outfmt", "11", "-out", args.output + ".asn"]
    print(" ".join(cline + args.blastn.split()))
    out = subprocess.call(cline + args.blastn.split())

    # output format xml (for parsing)
    cline = ["blast_formatter", "-archive", args.output + ".asn",
        "-outfmt", "5", "-out", args.output + ".xml"]
    print(" ".join(cline + args.blastn.split()))
    out = subprocess.call(cline)

    print("Completed in {:d} seconds.".format(int(round(time.time()-start))))

if __name__ == "__main__":
    main(sys.argv[1:])
