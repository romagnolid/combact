#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import os.path
import sys

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank",metavar="INPUT_FILE",
        help="the genbank input file")
    parser.add_argument("-o","--output",metavar="OUTPUT_FILE",
        help="the output fasta file [default stdout]")
    parser.add_argument("--igr", action="store_true",
        help="extract intergenic region from genbank file")

    args = parser.parse_args(args)

    if args.output and os.path.exists(args.output):
        raise IOError("{} already exists".format(args.output))
    elif args.output:
        output = open(args.output,"a")
    else:
        output = sys.stdout
    
    sequences = []
    records = SeqIO.parse(args.genbank,"gb")
    for record in records:
        start_igr = 0
        record_data = "gi|{gi}|ref|{ref}| {desc}".format(gi=record.annotations["gi"],ref=record.id,desc=record.description)
        for feat in record.features:
            if feat.type not in ("source","gene"):
                end_igr = feat.location.start
                if (end_igr - start_igr > 1) and args.igr:
                    seq_id = "igr|[{}:{}]".format(start_igr+1, end_igr)
                    seq = record.seq[start_igr:end_igr]
                    sequence = SeqRecord(seq,id=seq_id,description=record_data)
                    SeqIO.write(sequence, output, "fasta")
                tag = feat.qualifiers.get("locus_tag",["unknown_locus_tag"])[0]
                seq_id = "{}|[{}:{}]|{}".format(tag,feat.location.start+1,feat.location.end,feat.type)
                seq = feat.extract(record.seq)
                sequence = SeqRecord(seq,id=seq_id,description=record_data)
                SeqIO.write(sequence, output, "fasta")
                start_igr = feat.location.end
        end_igr = len(record.seq)
        if (end_igr - start_igr > 1) and args.igr:
            seq_id = "igr|[{}:{}]".format(start_igr+1, end_igr)
            seq = record.seq[start_igr:end_igr]
            sequence = SeqRecord(seq,id=seq_id,description=record_data)
            SeqIO.write(sequence, output, "fasta")

    if args.output:
        output.close()

if __name__ == "__main__":
    main(sys.argv[1:])
