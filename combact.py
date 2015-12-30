#!/usr/bin/env python
from __future__ import print_function, division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import argparse
import os
import os.path
import sys
import tools
import csv
import time

class AlignmentError(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        return repr(self.parameter)

def main(argv=None):
    parser = argparse.ArgumentParser(description="ComBact_0.6")
    parser.add_argument("input",metavar="INPUT_FILE",
        help="Blast output in xml format")
    parser.add_argument("-i","--inlist",metavar="GENOMES_LIST",default="input_list.txt",
        help="List of files [default=input_list.txt]")
    parser.add_argument("-o","--output", default="Output",
        help="output folder (default=Output)")
    parser.add_argument("-L","--length-cutoff", default=70, type=float,
        help="percent length cutoff [default=70]")
    parser.add_argument("-I","--identity-cutoff", default=80, type=float,
        help="percent identity cutoff [default=80]") 
    parser.add_argument("--silent-mut", action="store_true",
        help="additionally report silent mutations")
    
    args = parser.parse_args(argv)

    start = time.time()
    print("\nCompare Bacterial Genomes, current time:",time.strftime("%d/%m/%y at %H:%M:%S"))

    id_cutoff = args.identity_cutoff
    len_cutoff = args.identity_cutoff

    # output folder
    try:
        os.mkdir(args.output)
    except OSError:
        pass
    
    with open(args.inlist) as infile:
        reader = csv.reader(infile,delimiter="\t")
        genomes = [row[1] for row in reader]
    
    full_csv = open(os.path.join(args.output,"full.csv"),"w")
    full_writer = csv.DictWriter(full_csv, fieldnames=["Gene"] + genomes)
    full_writer.writeheader()
    
    nucl_csv = open(os.path.join(args.output,"nucl.csv"),"w")
    nucl_writer = csv.DictWriter(nucl_csv, fieldnames=["Gene"] + genomes)
    nucl_writer.writeheader()

    amino_csv = open(os.path.join(args.output,"amino.csv"),"w")
    amino_writer = csv.DictWriter(amino_csv, fieldnames=["Gene"] + genomes)        
    amino_writer.writeheader()
    
    blast_records = NCBIXML.parse(open(args.input))
    for blast_record in blast_records:
        query = blast_record.query
        q_len = blast_record.query_length
        #print(query)

        full_hits = {"Gene":query}
        nucl_hits = {"Gene":query}
        amino_hits = {"Gene":query}

        for alignment in blast_record.alignments:
            sbjct = alignment.hit_def
            
            for hsp in alignment.hsps:
                q_start = hsp.query_start
                q_end = hsp.query_end
                align_len = hsp.align_length
                identity = hsp.identities/align_len*100
                gaps = hsp.gaps
                q_seq = hsp.query
                s_seq = hsp.sbjct
     
                
                # wild-type
                if q_len == align_len and identity == 100:
                    #print("wt")
                    full_hits[sbjct] = full_hits.get(sbjct,[]) + ["wt"]
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + ["wt"]
                    amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["wt"]
    
                # SNP
                elif q_len == align_len and gaps == 0 and identity > id_cutoff:
                    #print("snp")
                    non_coding = tools.SNP_non_coding(q_seq,s_seq)
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + ["c.["+non_coding+"]"]
                    iscds = ("CDS" in query)

                    if iscds:
                        coding = tools.SNP_coding(q_seq,s_seq,True)
                        full_hits[sbjct] = full_hits.get(sbjct,[]) + ["c.["+non_coding+"]"+"="+"p.["+coding+"]"]
                        if args.silent_mut:
                            amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["p.["+coding+"]"]
                        else:
                            coding = tools.SNP_coding(q_seq,s_seq,False)
                            if len(coding)>0:
                                amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["p.["+coding+"]"]
                    else:
                        full_hits[sbjct] = full_hits.get(sbjct,[]) + ["c.["+non_coding+"]"]
                        amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["mutant"]
    
                # indel
                elif q_start == 1 and q_end == q_len and identity > id_cutoff: # and gaps > 0 (implicit)
                    #print("indel")
                    ins = tools.insertion(q_seq,s_seq)
                    dels = tools.deletion(q_seq,s_seq)

                    if ins and dels: # empty strings equal to false
                        indels = ins + ";" + dels
                    else:
                        indels = ins + dels

                    full_hits[sbjct] = full_hits.get(sbjct,[]) + ["c.["+indels+"]"]
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + ["c.["+indels+"]"]
                    amino_hits[sbjct] = amino_hits.get(sbjct,[]) + [indels]
    
                # fragment WT
                elif (q_start > 1 or q_end < q_len) and identity == 100 and align_len/q_len*100 > len_cutoff:
                    frag = "wt[{}:{}]".format(q_start,q_end)
                    full_hits[sbjct] = full_hits.get(sbjct,[]) + [frag]
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + [frag]
                    amino_hits[sbjct] = amino_hits.get(sbjct,[]) + [frag]

                # fragment WT
                elif q_start > 1 or q_end < q_len :
                    frag = "Fragment_mut[{}:{}]".format(q_start,q_end)
                    full_hits[sbjct] = full_hits.get(sbjct,[]) + [frag]
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + [frag]
                    amino_hits[sbjct] = amino_hits.get(sbjct,[]) + [frag]
    
                elif align_len/q_len*100 <= len_cutoff or identity <= id_cutoff:
                    #print("false_positive")
                    pass

                # any unforseen mutation
                else:
                    raise AlignmentError("Unknown mutation between " + query + " and " + sbjct)

        for genome in genomes:
            if full_hits.has_key(genome):
                full_hits[genome] = ";".join(full_hits[genome])
            if nucl_hits.has_key(genome):
                nucl_hits[genome] = ";".join(nucl_hits[genome])
            if amino_hits.has_key(genome):
                amino_hits[genome] = ";".join(amino_hits[genome])

        full_writer.writerow(full_hits)
        nucl_writer.writerow(nucl_hits)
        amino_writer.writerow(amino_hits)

    full_csv.close()
    nucl_csv.close()
    amino_csv.close()

    print("Completed in",round(time.time()-start,4),"seconds.")

if __name__ == "__main__":
    main(sys.argv[1:])
