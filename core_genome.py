#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
import argparse
import csv
import os
import os.path
import sys
import time

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("blastxml",metavar="INPUT_FILE",default="blast.xml",
        help="blast xml output to be parser")
    parser.add_argument("-o","--output",metavar="OUTPUT",default="Output",
        help="folder to store core genomes as fasta files [default=Output]")
    parser.add_argument("-i","--inlist",metavar="GENOMES_LIST",default="input_list.txt",
        help="List of files making the database [default=input_list.txt]")
    args = parser.parse_args(args)

    try:
        os.mkdir(args.output)
    except OSError:
        pass

    with open(args.inlist) as infile:
        reader = csv.reader(infile,delimiter="\t")
        genomes = [row[1] for row in reader]

    # overwrite files
    # this is needed because files will be opened with "append" option
    for genome in genomes:
        open(os.path.join(args.output,"{}_core.fa".format(genome)),"w+").read()

    # numeric table of core genomes
    csvfile = open("tag_table.csv","w")
    writer = csv.DictWriter(csvfile, fieldnames=["Locus tag"]+genomes)
    writer.writeheader()

    print("Parsing blast XML output, current time:",time.strftime("%d/%m/%y %H:%M:%S"),file=sys.stderr)
    start = time.time()
    # parse blastXML
    records = NCBIXML.parse(open(args.blastxml))
    # check whether a genome has a specific gene
    for record in records:
        hits = dict.fromkeys(genomes,0)
        for align in record.alignments:
            for hsp in align.hsps:
                if hsp.align_length >= record.query_length:
                    genome = align.hit_def.split()[-1] 
                    hits[genome] += 1
        # if present in all is core
        if all([value>0 for value in hits.values()]):
            #print(record.query,"is core",file=sys.stderr)
            identities = {}
            sequences = {}
            for align in record.alignments:
                genome = align.hit_def.split()[-1] 
                # keep only first hsp (best alignment)
                hsp = align.hsps[0]
                # keep only hsp with highest identity
                if hsp.identities > identities.get(genome,0):
                    identities[genome] = hsp.identities
                    # remove gaps
                    sbjct = hsp.sbjct.replace("-","") 
                    sequences[genome] = SeqRecord(Seq(sbjct),id=record.query)
            for genome,seq in sequences.items():
                with open(os.path.join(args.output,"{}_core.fa".format(genome)),"a") as core:
                    SeqIO.write(seq, core, "fasta")
        hits["Locus tag"] = record.query
        writer.writerow(hits)
    csvfile.close()
    print("Completed in",round(time.time()-start,4),"seconds.",file=sys.stderr)
    
if __name__ == "__main__":
    main(sys.argv[1:])
