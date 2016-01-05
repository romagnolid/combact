#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import os.path
import sys
import time

def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank",metavar="INPUT_FILE",
        help="the genbank input file")
    parser.add_argument("output",metavar="OUTPUT_FILE",
        help="the output fasta file")
    parser.add_argument("-t","--tag-list",metavar="FILE",
        help="list of tags to be extracted from genbank [optional]")        
    parser.add_argument("--igr", action="store_true",
        help="also extract intergenic regions from genbank [optional]")

    args = parser.parse_args(args)

    print("Parsing genbank, current time:",time.strftime("%d/%m/%y %H:%M:%S"),file=sys.stderr)
    start = time.time()

    if args.output:
        output = open(args.output,"w")
    else:
        output = sys.stdout
        
    # create multifasta with list of tag
    if args.tag_list:
        print("Proceed to extract tags from list",file=sys.stderr)
        tag_list = set(open(args.tag_list).read().split())
    else:
        print("Warning: tag-list not provided, all loci will be extracted",file=sys.stderr)
        
    if args.igr:
        print("Proceed to extract intergenic regions",file=sys.stderr)
    else:
        print("Warning: intergenic region will NOT be extracted",file=sys.stderr)
        
    sequences = []
    records = SeqIO.parse(args.genbank,"gb")
    for record in records:
        start_igr = 0
        record_data = "gi|{gi}|ref|{ref}| {desc}".format(gi=record.annotations["gi"],ref=record.id,desc=record.description)
        for feat in record.features:
            if feat.type not in ("source","gene") and feat.qualifiers.has_key("locus_tag"):
                end_igr = feat.location.start
                if (end_igr - start_igr > 1) and args.igr:
                    seq_id = "igr|[{s}:{e}]".format(s=start_igr+1, e=end_igr)
                    seq = record.seq[start_igr:end_igr]
                    sequence = SeqRecord(seq,id=seq_id,description=record_data)
                    SeqIO.write(sequence, output, "fasta")
                tag = feat.qualifiers["locus_tag"][0]
                seq_id = "{tag}|[{s}:{e}]|{typ}|".format(tag=tag,s=feat.location.start+1,e=feat.location.end,typ=feat.type)
                seq = feat.extract(record.seq)
                sequence = SeqRecord(seq,id=seq_id,description=record_data)

                if args.tag_list and (tag in tag_list):
                    SeqIO.write(sequence, output, "fasta")
                elif args.tag_list:
                    print("Warning:",tag,"[{typ}]".format(typ=feat.type),"not present",file=sys.stderr)
                else:
                    SeqIO.write(sequence, output, "fasta")                

                start_igr = feat.location.end
        
        end_igr = len(record.seq)       
        if (end_igr - start_igr > 1) and args.igr:
            seq_id = "igr|[{s}:{e}]".format(s=start_igr+1, e=end_igr)
            seq = record.seq[start_igr:end_igr]
            sequence = SeqRecord(seq,id=seq_id,description=record_data)
            SeqIO.write(sequence, output, "fasta")

    if args.output:
        output.close()
    print("Completed in",round(time.time()-start,4),"seconds.",file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])
