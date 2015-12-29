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


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta",metavar="FASTA",
        help="Multifasta of core genes")
    parser.add_argument("genomes",metavar="GENOMES",nargs="+",
        help="Files that will make the database")
    parser.add_argument("-o","--output",default="blast.xml",
        help="blast xml output [default blast.xml]")
    args = parser.parse_args(args)
    
    database = os.path.join(os.getcwd(),"Database","database.fa")
    try:
        os.mkdir("Database")
    except OSError:
        pass

    # collect all file paths
    filespaths = []
    for path in args.genomes:
        if os.path.isdir(path):
            filespaths += [os.path.join(path,x) for x in os.listdir(path)]
        elif os.path.isfile(path):
            filespaths += [path]
        else:
            raise IOError("{} is neither a folder or a file".format(path))

    with open("input_list.txt", "w") as glist:
        genomes = [os.path.splitext(os.path.basename(path))[0] for path in filespaths]
        glist.writelines(["{}\n".format(x) for x in genomes])
        
    # concatenate all fastas
    sequences = []
    for path in filespaths:
        genomename = os.path.splitext(os.path.basename(path))[0]
        records = SeqIO.parse(path,"fasta")
        for i,record in enumerate(records,1):
            myID = "{} ({})".format(record.id,i)
            sequences.append(SeqRecord(record.seq,id=myID,
                description=genomename))
    SeqIO.write(sequences, database,"fasta")

    subprocess.call(["makeblastdb","-parse_seqids","-dbtype","nucl",
        "-in", database,"-out", os.path.splitext(database)[0]])

    # blast multifasta
    cline = NcbiblastnCommandline(query=args.fasta, db=os.path.splitext(database)[0], outfmt=5, out=args.output)
    print("\n",time.asctime(),sep="")
    print(str(cline))
    cline()
    print(time.asctime())        
            
if __name__ == "__main__":
    main(sys.argv[1:])
