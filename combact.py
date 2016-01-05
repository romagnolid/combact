#!/usr/bin/env python
from __future__ import print_function, division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import argparse
import os
import os.path
import sys
import csv
import time

def SNP_coding(x,y,report_silent=False):
    """Return string differences as SNP mutation at protein level"""
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    start_table = {"TTG":"M","CTG":"M","ATT":"M","ATC":"M","ATA":"M","ATG":"M","GTG":"M"}

    mutations = []
    x = x.upper()
    y = y.upper()
    x_codons = [x[i:i+3] for i in range(0, len(x)-3+1, 3)]
    y_codons = [y[i:i+3] for i in range(0, len(y)-3+1, 3)]

    for i in range(len(x_codons)):
        if x_codons[i] != y_codons[i]:
            if i == 0:
                a = start_table.get(x_codons[i],"X")
                b = start_table.get(y_codons[i],"X")
            else:
                a = codon_table.get(x_codons[i],"X")
                b = codon_table.get(y_codons[i],"X")

            if a == b != "X" and report_silent:
                mutations.append("{}{}{}".format(a, i+1, b))
            elif a != b == "*":
                mutations.append("{}{}{}".format(a, i+1, b))
            elif a != b or a == b == "X":
                mutations.append("{}{}{}".format(a, i+1, b))
    return(";".join(mutations))

def SNP_non_coding(x,y):
    """Return string differences as SNP mutation at DNA level"""
    mutations = []
    x = x.upper()
    y = y.upper()
    i = 0
    while i < len(x):
        if x[i] != y[i]:
            mutated = True
            j = i
            i += 1
            while mutated and i<len(x):
                if x[i] != y[i]:
                    i += 1
                else:
                    mutations.append("{}{}>{}".format(j+1,x[j:i],y[j:i]))
                    mutated = False
                    i += 1
            if i == len(x):
                mutations.append("{}{}>{}".format(j+1,x[j:i],y[j:i]))
        else:
            i += 1
    return(";".join(mutations))

def insertion(x,y):
    mutations = []
    x = x.upper()
    y = y.upper()
    i = 0
    l = 0 # nucleotides count (excluding gaps)
    while i < len(x):
        if x[i] == "-":
            j = i
            has_gap = True
            i += 1
            while has_gap:
                if x[i] != "-":
                    mutations.append("{}_{}ins{}".format(l,l+1, y[j:i]))
                    has_gap = False
                    i += 1
                    l += 1
                else:
                    i += 1
        else:
            i += 1
            l += 1
    return(";".join(mutations))

def deletion(x,y):
    mutations = []
    x = x.upper()
    y = y.upper()
    i = 0
    k = 0
    while i < len(y):
        if y[i] == "-":
            j = i
            has_gap = True
            i += 1
            while has_gap:
                if y[i] != "-":
                    if (i-j) > 1:
                        mutations.append(
                            "{}_{}del{}".format(j+1-k,j+(i-j)-k,x[j:i]))
                    else:
                        mutations.append(
                            "{}del{}".format(j+1-k,x[j:i]))
                    if x[i] == "-":
                        k += 1
                    has_gap = False
                    i += 1
                else:
                    i += 1
        else:
            if x[i] == "-":
                k += 1
            i += 1
    return(";".join(mutations))

def main(argv=None):
    parser = argparse.ArgumentParser(description="ComBact_0.6.2")
    parser.add_argument("input",metavar="INPUT_FILE",
        help="the blast xml output")
    parser.add_argument("output",metavar="OUTPUT_DIRECTORY",
        help="the output directory containing the three csv files")
    parser.add_argument("-i","--inlist",metavar="GENOMES_LIST",default="input_list.txt",
        help="list of genomes used to make the database [default=\"input_list.txt\"]")
    parser.add_argument("--silent-mut", action="store_true",
        help="additionally report silent mutations")
    parser.add_argument("-L","--length",metavar="CUTOFF",default=70,type=float,dest="len_cutoff",
        help="percent alignment length cutoff [default=70]")
    parser.add_argument("-I","--identity",metavar="CUTOFF",default=80,type=float,dest="id_cutoff",
        help="percent identity cutoff [default=80]")

    args = parser.parse_args(argv)

    start = time.time()
    print("Compare Bacterial Genomes, current time:",time.strftime("%d/%m/%y %H:%M:%S"))

    genomes = list()
    # list of files used for the database
    with open(args.inlist) as infile:
        reader = csv.reader(infile,delimiter="\t")
        for row in reader:
            if len(row) == 1:
                genomes.append(os.path.basename(os.path.splitext(row[0])[0]))
            else:
                genomes.append(row[1])

    # output folder
    try:
        os.mkdir(args.output)
    except OSError:
        pass

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
        full_hits = {"Gene":query}
        nucl_hits = {"Gene":query}
        amino_hits = {"Gene":query}
        for alignment in blast_record.alignments:
            sbjct = alignment.hit_def.split()[-1]
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
                elif q_len == align_len and gaps == 0 and identity > args.id_cutoff:
                    #print("snp")
                    non_coding = SNP_non_coding(q_seq,s_seq)
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + ["c.["+non_coding+"]"]
                    iscds = ("CDS" in query)
                    if iscds:
                        coding = SNP_coding(q_seq,s_seq,True)
                        full_hits[sbjct] = full_hits.get(sbjct,[]) + ["c.["+non_coding+"]"+"="+"p.["+coding+"]"]
                        if args.silent_mut is True:
                            amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["p.["+coding+"]"]
                        else:
                            coding = SNP_coding(q_seq,s_seq,False)
                            if len(coding)>0:
                                amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["p.["+coding+"]"]
                    else:
                        full_hits[sbjct] = full_hits.get(sbjct,[]) + ["c.["+non_coding+"]"]
                        amino_hits[sbjct] = amino_hits.get(sbjct,[]) + ["mutant"]
                # indel
                elif q_start == 1 and q_end == q_len and identity > args.id_cutoff: # and gaps > 0 (implicit)
                    #print("indel")
                    ins = insertion(q_seq,s_seq)
                    dels = deletion(q_seq,s_seq)

                    if ins and dels: # empty strings equal to false
                        indels = ins + ";" + dels
                    else:
                        indels = ins + dels

                    full_hits[sbjct] = full_hits.get(sbjct,[]) + ["c.["+indels+"]"]
                    nucl_hits[sbjct] = nucl_hits.get(sbjct,[]) + ["c.["+indels+"]"]
                    amino_hits[sbjct] = amino_hits.get(sbjct,[]) + [indels]

                # fragment WT
                elif (q_start > 1 or q_end < q_len) and identity == 100 and align_len/q_len*100 > args.len_cutoff:
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

                elif align_len/q_len*100 <= args.len_cutoff or identity <= args.id_cutoff:
                    #print("false_positive")
                    pass

                # any unforseen mutation
                else:
                    print("Error: unknown mutation between",query,"and",sbjct,file=sys.stderr)
                    sys.exit()

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
