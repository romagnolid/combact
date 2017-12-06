#!/usr/bin/python

from __future__ import print_function, division
import argparse
import os
import os.path
import sys
import csv
import time
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML

from itertools import groupby
from operator import itemgetter

def snp_coding(x, y, k=1, print_silent=False):
    """Return a string of SNPs at protein level"""

    mutations = []
    x_amino = Seq(x, generic_dna).translate(table=11)
    y_amino = Seq(y, generic_dna).translate(table=11)
    # do not use cds option: it removes stop codon *;
    # it throws error both if first codon is not a start codon or last codon is not end codon
    # better check if first codon is start codon ecc.

    x_codons = [x[i:i+3] for i in range(0, len(x)-3+1, 3)]
    y_codons = [y[i:i+3] for i in range(0, len(y)-3+1, 3)]

    indexes = [i for i in range(len(x_codons)) if x_codons[i] != y_codons[i]]
    for h, g in groupby(enumerate(indexes), lambda (i, x): i-x):
        ind = map(itemgetter(1), g)

        # silent
        for i in ind:
            if x_amino[i] == y_amino[i] and print_silent:
                mutations.append("{x}{i}{y}".format(i=i+k, x=x_amino[i], y=y_amino[i]))
            else:
                mutations.append("{x}{i}{y}".format(i=i+k, x=x_amino[i], y=y_amino[i]))

    return(";".join(mutations))

def snp_non_coding(x, y, k=1):
    """Return a string of SNPs at DNA level.
    x = seq1
    y = seq2
    k = starting nucleotide default = 1"""

    mutations = []
    x = x.upper()
    y = y.upper()
    indexes = [i for i in range(len(x)) if x[i] != y[i]]
    for (h, g) in groupby(enumerate(indexes), lambda (i, x): i-x):
        i = map(itemgetter(1), g)
        if len(i) == 1:
            # e.g. "66T>G"
            mutations.append(
                "{i}{x}>{y}".format(i=i[0]+k, x=x[i[0]], y=y[i[0]]))
        else:
            # e.g. "c.112_117delAGGTCAinsTG"
            mutations.append(
                "{i}_{j}del{x}ins{y}".format(i=i[0]+k, j=i[-1]+k, x=x[i[0]:i[-1]+1], y=y[i[0]:i[-1]+1]))
    return(";".join(mutations))

def insertion(x, y):
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
                    # e.g. "c.112_113insAGA"
                    mutations.append("{}_{}ins{}".format(l, l+1, y[j:i]))
                    has_gap = False
                    i += 1
                    l += 1
                else:
                    i += 1
        else:
            i += 1
            l += 1
    return(";".join(mutations))

def deletion(x, y):
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
                            "{}_{}del{}".format(j+1-k, j+(i-j)-k, x[j:i]))
                    else:
                        mutations.append(
                            "{}del{}".format(j+1-k, x[j:i]))
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
    parser = argparse.ArgumentParser(description="ComBact-0.8")
    parser.add_argument("input", metavar="INPUT_FILE",
        help="XML blast alignment")
    parser.add_argument("output", metavar="OUTPUT_DIRECTORY",
        help="output directory containing three tsv files")
    parser.add_argument("-i", "--input_list", metavar="GENOMES_LIST",
        help="list of genome alignments to be parsed")
    parser.add_argument("--silent_mut", action="store_true",
        help="report silent mutations")
    parser.add_argument("-L", "--length", metavar="CUTOFF", default=0, type=float, dest="len_cutoff",
        help="percent alignment length cutoff [default=0]")
    parser.add_argument("-I", "--identity", metavar="CUTOFF", default=0, type=float, dest="id_cutoff",
        help="percent identity cutoff [default=0]")

    args = parser.parse_args(argv)

    start = time.time()
    print("Filters:\n  identity = {id}\n  length = {len}".format(
        id=args.id_cutoff, len=args.len_cutoff))
    print("Reporting silent mutations: {}".format(args.silent_mut))

    print("\nComBact, current time:", time.strftime("%d/%m/%y %H:%M:%S"))
    with open(args.input_list) as infile:
        genomes = [os.path.basename(os.path.splitext(line.split()[0])[0])
            for line in infile]

    # output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    full_tsv = open(os.path.join(args.output, "full.tsv"), "w")
    full_writer = csv.DictWriter(full_tsv, delimiter="\t",
        fieldnames=["Gene"] + genomes)
    full_writer.writeheader()

    nucl_tsv = open(os.path.join(args.output, "nucl.tsv"), "w")
    nucl_writer = csv.DictWriter(nucl_tsv, delimiter="\t",
        fieldnames=["Gene"] + genomes)
    nucl_writer.writeheader()

    amino_tsv = open(os.path.join(args.output, "amino.tsv"), "w")
    amino_writer = csv.DictWriter(amino_tsv, delimiter="\t",
        fieldnames=["Gene"] + genomes)
    amino_writer.writeheader()

    blast_records = NCBIXML.parse(open(args.input))
    # one blast_record for each sequence in reference
    for blast_record in blast_records:
        query = blast_record.query
        q_len = blast_record.query_length

        full_hits  = {"Gene": query}
        nucl_hits  = {"Gene": query}
        amino_hits = {"Gene": query}

        # one alignment for each different fasta sequence in subjects:
        # genomes may result in multiple alignments if splitted in different fasta sequences
        for alignment in blast_record.alignments:
            sbjct = alignment.accession
            if (sbjct in genomes):
            # if (sbjct not in nucl_hits):
                hsp = alignment.hsps[0] #high scoring pair
                q_start = hsp.query_start
                q_end = hsp.query_end
                align_len = hsp.align_length
                identity = hsp.identities/align_len*100
                gaps = hsp.gaps
                q_seq = hsp.query
                s_seq = hsp.sbjct

                # wild-type
                if q_len == align_len and identity == 100:
                    full_hits[sbjct]  = "wt"
                    nucl_hits[sbjct]  = "wt"
                    amino_hits[sbjct] = "wt"

                # SNP
                elif (q_len == align_len and gaps == 0 and identity > args.id_cutoff):
                    non_coding = snp_non_coding(q_seq, s_seq)
                    nucl_hits[sbjct] = "c.[{}]".format(non_coding)

                    if ("CDS" in query):
                        if (len(q_seq) % 3 != 0) or (len(s_seq) % 3 != 0):
                            print(query, "length is not multiple of three")
                        coding = snp_coding(q_seq, s_seq, True)
                        full_hits[sbjct] = "c.[{}]=p.[{}]".format(non_coding, coding)

                        if (args.silent_mut):
                            amino_hits[sbjct] = "p.[{}]".format(coding)
                        else:
                            coding = snp_coding(q_seq, s_seq, False)
                            if (len(coding) > 0):
                                amino_hits[sbjct] = "p.[{}]".format(coding)
                    else:
                        full_hits[sbjct] = "c.[{}]".format(non_coding)
                        amino_hits[sbjct] = "mutant"

                # indel
                elif (q_start == 1 and q_end == q_len and identity > args.id_cutoff): # and gaps > 0 (implicit)
                    ins = insertion(q_seq, s_seq)
                    dels = deletion(q_seq, s_seq)

                    if ins and dels: # empty strings equal to false
                        indels = ins + ";" + dels
                    else:
                        indels = ins + dels

                    full_hits[sbjct]  = "c.[{}]".format(indels)
                    nucl_hits[sbjct]  = "c.[{}]".format(indels)
                    amino_hits[sbjct] = "indels"

                # fragment WT
                elif (q_start > 1 or q_end < q_len) and identity == 100 and align_len/q_len*100 > args.len_cutoff:
                    frag = "wt[{}:{}]".format(q_start, q_end)
                    full_hits[sbjct]  = frag
                    nucl_hits[sbjct]  = frag
                    amino_hits[sbjct] = frag

                # fragment mutant
                elif (q_start > 1 or q_end < q_len) and align_len/q_len*100 > args.len_cutoff:
                    # snps only
                    if gaps == 0:
                        #print(query)
                        non_coding = snp_non_coding(q_seq, s_seq, k=q_start)
                        nucl_hits[sbjct] = "c.({}:{})[{}]".format(q_start, q_end, non_coding)
                        full_hits[sbjct] = "c.({}:{})[{}]".format(q_start, q_end, non_coding)

                elif align_len/q_len*100 <= args.len_cutoff or identity <= args.id_cutoff:
                    pass

                else: # debugging purpose: any unforseen mutation
                    print("Error: unknown mutation between", query, "and", sbjct, file=sys.stderr)
                    sys.exit()

        full_writer.writerow(full_hits)
        nucl_writer.writerow(nucl_hits)
        amino_writer.writerow(amino_hits)

    full_tsv.close()
    nucl_tsv.close()
    amino_tsv.close()

    print("Completed in {:d} seconds.".format(int(round(time.time()-start))))

if __name__ == "__main__":
    main(sys.argv[1:])
