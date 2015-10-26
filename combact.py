#!/usr/bin/python
from __future__ import print_function, division
from subprocess import call
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from sys import exit
import os.path
import tools

def print_fasta(sequence,handle=None):
    i = 0
    while i < len(sequence):
        print(sequence[i:i+80], file=handle)
        i += 80

def write_igr(locations, record, handle=None):
    """Write intergenic regions into an open file
    """
    if locations[0].strand == 1:
        strand = "+"
    else:
        strand = "-"
    i = 0
    start = 0
    while i < len(locations):
        end = locations[i].start
        if end - start > 0:
            header = ">IGR:[{}:{}]({})|{}".format(start, end, strand, record.name)
            seq = record.seq[start:end]
            print(header,file=handle)
            print_fasta(seq, handle)
        start = locations[i].end
        i += 1
    end = len(record.seq)
    if end - start > 0:
        header = ">IGR:[{}:{}]({})|{}".format(start, end, strand, record.name)
        seq = record.seq[start:end]
        print(header,file=handle)
        print_fasta(seq, handle)
            
def parse_gb_record(infilename,handle=None):
    """Read genebank file and convert it to a fasta file
    """
    gb_records = SeqIO.parse(infilename,"gb")
    genes_location_plus = []
    genes_location_minus = []

    for gb_record in gb_records:
        for feature in gb_record.features:
            if feature.type not in ("gene","source","variation"):
                locus_tag = feature.qualifiers.get("locus_tag")
                if locus_tag:
                    header = ">{}:{}:{}".format(feature.type,feature.location,locus_tag[0])
                else:
                    header = ">{}:{}".format(feature.type,feature.location)
                seq = feature.extract(gb_record.seq)
                print(header,file=handle)
                print_fasta(seq, handle)
                #db_xref = feature.qualifiers.get("db_xref",[])
                #gene = feature.qualifiers.get("gene",[])
                #product = feature.qualifiers.get("product",[])
            elif feature.type == "gene" and feature.strand == 1:
                genes_location_plus.append(feature.location)
            elif feature.type == "gene" and feature.strand == -1:
                genes_location_minus.append(feature.location)
        
        write_igr(genes_location_plus, gb_record, handle)
        write_igr(genes_location_minus, gb_record, handle)

def concatenate_fasta(infilenames, handle=None):
    """ Concatenate several multifasta into a single multifasta
    """
    if not hasattr(infilenames,"__iter__"):
        raise TypeError("First argument must be iterable")

    for path in infilenames:
        filename = os.path.splitext(os.path.basename(path))[0]
        for seq_record in SeqIO.parse(path,"fasta"):
            header = ">{}:{}".format(filename, seq_record.id)
            print(header,file=handle)
            print_fasta(seq_record.seq, handle)

def makeblastdb(db_file):
    call(["makeblastdb","-in", db_file,"-dbtype","nucl","-out", db_file.split(".")[0]])

def blastn(query, db_file, out):
    cline = NcbiblastnCommandline(query=query, db=db_file.split(".")[0], outfmt=5, out=out)
    cline()

def get_mutations(infilename, inputlist, handle_nucl=None, handle_amino=None, length_cutoff=70, identity_cutoff=80, with_silent=False):
    sep="\t"
    header = sep.join(["GeneID"]+inputlist)
    print(header,file=handle_nucl)
    print(header,file=handle_amino)

    blast_records = NCBIXML.parse(open(infilename))

    for blast_record in blast_records:
        query = blast_record.query
        q_len = int(blast_record.query_length)
        nucl_array = dict(zip(inputlist, [[] for i in range(len(inputlist))]))
        amino_array = dict(zip(inputlist, [[] for i in range(len(inputlist))]))
        for alignment in blast_record.alignments:
            sbjct = alignment.hit_def.split(":")[0]
    
            for hsp in alignment.hsps:
                q_start = int(hsp.query_start)
                q_end = int(hsp.query_end)
                align_len = int(hsp.align_length)
                identity = int(hsp.identities)/align_len*100
                gaps = int(hsp.gaps)
                q_seq = hsp.query
                s_seq = hsp.sbjct
     
                if identity < identity_cutoff or align_len/q_len*100 < length_cutoff:
                    nucl_array[sbjct].append("non_functional")
                
                # wild-type
                elif identity == 100 and q_len == align_len:
                    nucl_array[sbjct].append("WT")
    
                # SNP
                elif identity < 100 and q_len == align_len and gaps == 0:
                    non_coding = tools.SNP_non_coding(q_seq,s_seq)
                    if with_silent:
                        coding = tools.SNP_coding(q_seq,s_seq,True)
                    else: 
                        coding = tools.SNP_coding(q_seq,s_seq)

                    # whether the gene codes for a protein or not
                    if query.startswith("CDS") or "peg" in query:
                        nucl_array[sbjct].append(non_coding+"$"+coding)
                        amino_array[sbjct].append(coding)
                    else:
                        nucl_array[sbjct].append(non_coding)
    
                # indel
                elif q_start == 1 and q_end == q_len and gaps > 0:
                    ins = tools.insertion(q_seq,s_seq)
                    dels = tools.deletion(q_seq,s_seq)
                    if len(ins)*len(dels) > 0: 
                        indels = ins+";"+dels
                    else:
                        indels = ins + dels
                    nucl_array[sbjct].append(indels)
                    amino_array[sbjct].append("INDEL")
    
                # fragment
                elif q_start > 1 or q_end < q_len:
                    frag = "Fragment[{0}:{1}]".format(q_start,q_end)
                    nucl_array[sbjct].append(frag)
    
                # any unforseen mutation
                else:
                    exit("Unknown_mutation")

        print(sep.join([query]+["/".join(nucl_array[sbjct]) for sbjct in inputlist]),file=handle_nucl)

        if sum([len(amino_array[sbjct]) for sbjct in inputlist]) > 0:
            print(sep.join([query]+["/".join(amino_array[sbjct]) for sbjct in inputlist]),file=handle_amino)
