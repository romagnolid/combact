#!/usr/bin/env python
from __future__ import print_function,division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from os.path import basename, splitext
from sys import stdout
import tools

class AlignmentError(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        return repr(self.parameter)

def parse_gb_file(infilename, handle_out=stdout, add_igr=False):
    """Read genebank file and convert it into a multifasta.
    """
    record = SeqIO.read(infilename,"gb")
    sequences = []

    start_igr = 0
    for feature in record.features:
        if feature.type in ("CDS","tRNA","rRNA"):
               end_igr = feature.location.start
               if end_igr - start_igr > 1 and add_igr:
                   seq_id = "IGR:[{}:{}]".format(start_igr + 1, end_igr)
                   seq = record.seq[start_igr:end_igr]
                   sequences.append(SeqRecord(seq, id=seq_id, description="intergenic region"))

               locus_tag = feature.qualifiers["locus_tag"][0]
               seq_desc = feature.qualifiers["product"][0]
               seq_id = "{}:[{}:{}]:{}".format(feature.type,feature.location.start + 1,feature.location.end,locus_tag)
               seq = feature.extract(record.seq)
               sequences.append(SeqRecord(seq,id=seq_id,description=seq_desc))
               start_igr = feature.location.end

    end_igr = len(record.seq)
    if end_igr - start_igr > 1 and add_igr:
        seq_id = "IGR:[{}:{}]".format(start_igr + 1, end_igr)
        seq = record.seq[start_igr:end_igr]
        sequences.append(SeqRecord(seq, id=seq_id, description="intergenic region"))
    SeqIO.write(sequences, handle_out, "fasta")

def concatenate_fasta(infilenames, handle_out=stdout):
    """Concatenate several fasta files into a single multifasta.
    """
    if not hasattr(infilenames,"__iter__"):
        raise TypeError("First argument must be iterable")

    sequences = []
    for path in infilenames:
        filename = splitext(basename(path))[0]
        records = SeqIO.parse(path,"fasta")
        for i, seq_record in enumerate(records, 1):
            seq_id = "{}|{}|{}".format(filename, seq_record.id, str(i)) # remove str(i)?
            sequences.append(SeqRecord(seq_record.seq, id=seq_id))
    SeqIO.write(sequences, handle_out ,"fasta")

def get_mutations(infilename, inputlist, handle_full=None, handle_nucl=None, handle_amino=None, identity_cutoff=80, length_cutoff=70, add_silent=False):
    sep="\t"
    header = sep.join(["GeneName"] + inputlist)
    print(header,file=handle_full)
    print(header,file=handle_nucl)
    print(header,file=handle_amino)

    blast_records = NCBIXML.parse(open(infilename))

    for blast_record in blast_records:
        query = blast_record.query
        q_len = blast_record.query_length
        full_array = dict(zip(inputlist, [[] for i in range(len(inputlist))]))
        nucl_array = dict(zip(inputlist, [[] for i in range(len(inputlist))]))
        amino_array = dict(zip(inputlist, [[] for i in range(len(inputlist))]))
        for alignment in blast_record.alignments:
            sbjct = alignment.hit_def.split("|")[0]
    
            for hsp in alignment.hsps:
                q_start = hsp.query_start
                q_end = hsp.query_end
                align_len = hsp.align_length
                identity = hsp.identities/align_len*100
                gaps = hsp.gaps
                q_seq = hsp.query
                s_seq = hsp.sbjct
     
                # short segment: false positive
                if identity <= identity_cutoff or align_len/q_len*100 <= length_cutoff:
                    pass
                
                # wild-type
                elif identity == 100 and q_len == align_len:
                    nucl_array[sbjct].append("WT")
    
                # SNP
                elif identity < 100 and gaps == 0:
                    non_coding = tools.SNP_non_coding(q_seq,s_seq)
                    if add_silent:
                        coding = tools.SNP_coding(q_seq,s_seq,True)
                    else: 
                        coding = tools.SNP_coding(q_seq,s_seq)

                    # whether the gene codes for a protein or not
                    full_array[sbjct].append(non_coding+"$"+coding)
                    nucl_array[sbjct].append(non_coding)
                    amino_array[sbjct].append(coding)
    
                # indel
                elif q_start == 1 and q_end == q_len and gaps > 0:
                    ins = tools.insertion(q_seq,s_seq)
                    dels = tools.deletion(q_seq,s_seq)
                    if len(ins)*len(dels) > 0: 
                        indels = ins + ";" + dels
                    else:
                        indels = ins + dels
                    full_array[sbjct].append(indels)
                    nucl_array[sbjct].append(indels)
                    amino_array[sbjct].append("INDEL")
    
                # fragment
                elif q_start > 1 or q_end < q_len:
                    frag = "Fragment[{0}:{1}]".format(q_start,q_end)
                    nucl_array[sbjct].append(frag)
    
                # any unforseen mutation
                else:
                    raise AlignmentError("Unknown mutation type occuring between " + query + " and " + sbjct)

        print(sep.join([query]+["/".join(full_array[sbjct]) for sbjct in inputlist]),file=handle_full)
        print(sep.join([query]+["/".join(nucl_array[sbjct]) for sbjct in inputlist]),file=handle_nucl)

        if sum([len(amino_array[sbjct]) for sbjct in inputlist]) > 0:
            print(sep.join([query]+["/".join(amino_array[sbjct]) for sbjct in inputlist]),file=handle_amino)
