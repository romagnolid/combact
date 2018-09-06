#!/usr/bin/env python2

from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import os.path
import sys
import time

strand_dict = {1:"+", -1:"-"}

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", metavar="INPUT_FILE",
        help="the genbank input file")
    parser.add_argument("output", metavar="OUTPUT_FILE",
        help="the output fasta file")
    parser.add_argument("-l", "--list", metavar="FILE",
        help="parse a selection of genes with locus tag list")
    parser.add_argument("--igr", action="store_true",
        help="extract intergenic regions")
    parser.add_argument("--embl", action="store_true",
        help="parse from embl flat format file instead of genbank")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.8")

    args = parser.parse_args(args)

    print("Parsing genbank/embl input file, current time:", time.strftime("%d/%m/%y %H:%M:%S"))
    start = time.time()

    if args.output:
        output = open(args.output, "w")
    else:
        output = sys.stdout

    if args.list:
        print("Extract tags from list")
        tag_list = set(open(args.list).read().split())
    else:
        print("Extract all loci")

    print("Extract intergenic regions: {}".format(args.igr))

    sequences = []
    if args.embl:
        records = SeqIO.parse(args.input, "embl")
    else:
        records = SeqIO.parse(args.input, "gb")

    for record in records:
        start_igr = 0
        n_igr = 1
        for feat in record.features:
            if feat.type in ("CDS", "tRNA", "rRNA"):
                # IGR
                end_igr = feat.location.start
                if (end_igr - start_igr > 0) and args.igr:
                    seq_id = "igr_{n:0>5}|{s}:{e}".format(n=n_igr, s=start_igr+1, e=end_igr)
                    seq = record.seq[start_igr:end_igr]
                    sequence = SeqRecord(seq, id=seq_id, description="")
                    SeqIO.write(sequence, output, "fasta")
                    n_igr += 1

                # GENE
                seq_id = "{tag}|{s}:{e}:{strand}|{typ}| {gene_name}:".format(
                    strand=strand_dict[feat.strand],
                    tag=feat.qualifiers.get("locus_tag", ["unknown_tag"])[0],
                    s=feat.location.start+1,
                    e=feat.location.end,
                    typ=feat.type,
                    gene_name=feat.qualifiers.get("gene", ["unknown_gene"])[0])
                seq = record.seq[feat.location.start:feat.location.end]

                sequence = SeqRecord(seq, id=seq_id,
                    description=feat.qualifiers.get("product", ["unknown_product"])[0])

                if args.list and (tag in tag_list):
                    SeqIO.write(sequence, output, "fasta")
                elif not args.list:
                    SeqIO.write(sequence, output, "fasta")

                start_igr = feat.location.end

        end_igr = len(record.seq)
        if (end_igr - start_igr > 0) and args.igr:
            seq_id = "igr_{n:0>5}|{s}:{e}".format(n=n_igr, s=start_igr+1, e=end_igr)
            seq = record.seq[start_igr:end_igr]
            sequence = SeqRecord(seq, id=seq_id, description="")
            SeqIO.write(sequence, output, "fasta")

    if args.output:
        output.close()
    print("Completed in {:d} seconds.".format(int(round(time.time()-start))))

if __name__ == "__main__":
    main(sys.argv[1:])
