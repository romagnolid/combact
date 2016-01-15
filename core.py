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
    parser = argparse.ArgumentParser(description="Core_0.1.0")
    parser.add_argument("blastxml",metavar="INPUT_FILE",default="blast.xml",
        help="blast xml output to be parser")
    parser.add_argument("output",metavar="OUTPUT_DIRECTORY",
        help="the output directory containing fasta files of core genomes")
    parser.add_argument("-i","--inlist",metavar="GENOMES_LIST",default="inlist.txt",
        help="List of files making the database [default=inlist.txt]")
    args = parser.parse_args(args)

    try:
        os.mkdir(args.output)
    except OSError:
        for path in os.listdir(args.output):
            os.remove(os.path.join(args.output,path))

    # list of files used for the database
    with open(args.inlist) as infile:
        genomes = [os.path.basename(os.path.splitext(line.split()[0])[0]) for line in infile]

    # numeric table of core genomes
    core_genome_csv = open("tag_table.csv","w")
    writer = csv.DictWriter(core_genome_csv, fieldnames=["Locus tag"]+genomes)
    writer.writeheader()

    print("Parsing blast XML output, current time:",time.strftime("%d/%m/%y %H:%M:%S"),file=sys.stderr)
    start = time.time()
    # parse blastXML
    records = NCBIXML.parse(open(args.blastxml))
    # check whether a genome has a specific gene
    for record in records:
        #print("Parsing {rec} record".format(rec=record.query),file=sys.stderr)
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
                with open(os.path.join(args.output,"{}.core.fa".format(genome)),"a") as core:
                    SeqIO.write(seq, core, "fasta")
        hits["Locus tag"] = record.query
        writer.writerow(hits)
    core_genome_csv.close()
    print("Completed in",round(time.time()-start,4),"seconds.",file=sys.stderr)
    
if __name__ == "__main__":
    main(sys.argv[1:])
