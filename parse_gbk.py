#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import os.path
import sys
import time

strand_dict = {1:"+",-1:"-"}

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank",metavar="INPUT_FILE",
        help="the genbank input file")
    parser.add_argument("output",metavar="OUTPUT_FILE",
        help="the output fasta file")
    parser.add_argument("-l","--list",metavar="FILE",
        help="selective parsing of genes with locus tag list")        
    parser.add_argument("--igr", action="store_true",
        help="extract intergenic regions")
    parser.add_argument("--embl", action="store_true",
        help="parse from embl flat format file instead of genbank")

    args = parser.parse_args(args)

    print("Parsing genbank, current time:",time.strftime("%d/%m/%y %H:%M:%S"))
    start = time.time()

    if args.output:
        output = open(args.output,"w")
    else:
        output = sys.stdout
        
    # create multifasta with list of tag
    if args.list:
        print("Proceed to extract tags from list")
        tag_list = set(open(args.list).read().split())
    else:
        print("Warning: tag-list not provided, all loci will be extracted")
        
    if args.igr:
        print("Proceed to extract intergenic regions")
    else:
        print("Warning: intergenic region will NOT be extracted")
        
    sequences = []
    if args.embl:
        records = SeqIO.parse(args.genbank,"embl")
    else:
        records = SeqIO.parse(args.genbank,"gb")

    for record in records:
        start_igr = 0
        for feat in record.features:
            if feat.type in ("CDS","tRNA","rRNA"):# and feat.qualifiers.has_key("locus_tag"):
                # IGR
                end_igr = feat.location.start
                if (end_igr - start_igr > 0) and args.igr:
                    seq_id = "igr|{s}:{e}".format(s=start_igr+1, e=end_igr)
                    #pre = record.seq[start_igr-9:start_igr].lower()
                    #post = record.seq[end_igr:end_igr+9].lower()
                    #seq = pre + record.seq[start_igr:end_igr] + post
                    seq = record.seq[start_igr:end_igr]
                    sequence = SeqRecord(seq,id=seq_id,description="")
                    SeqIO.write(sequence, output, "fasta")

                # GENE
                seq_id = "{tag}|{s}:{e}:{strand}|{typ}| {gene_name}:".format(
                    strand=strand_dict[feat.strand], 
                    tag=feat.qualifiers.get("locus_tag",["unknown_tag"])[0], 
                    s=feat.location.start+1, 
                    e=feat.location.end, 
                    typ=feat.type,
                    gene_name=feat.qualifiers.get("gene",["unknown_gene"])[0])
                #pre = record.seq[feat.location.start-9:feat.location.start].lower()
                #post = record.seq[feat.location.end:feat.location.end+9].lower()
                seq = record.seq[feat.location.start:feat.location.end]
                #seq = feat.extract(record.seq)

                sequence = SeqRecord(seq, id=seq_id,
                    description=feat.qualifiers.get("product",["unknown_product"])[0])

                if args.list and (tag in tag_list):
                    SeqIO.write(sequence, output, "fasta")
                elif args.list and not (tag in tag_list):
                    pass
                elif not args.list:
                    SeqIO.write(sequence, output, "fasta")                

                start_igr = feat.location.end
        
        end_igr = len(record.seq)       
        if (end_igr - start_igr > 0) and args.igr:
            seq_id = "igr|{s}:{e}".format(s=start_igr+1, e=end_igr)
            seq = record.seq[start_igr:end_igr]
            sequence = SeqRecord(seq,id=seq_id,description="")
            SeqIO.write(sequence, output, "fasta")

    if args.output:
        output.close()
    print("Completed in",round(time.time()-start,4),"seconds.")

if __name__ == "__main__":
    main(sys.argv[1:])
