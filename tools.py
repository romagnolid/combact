#!/usr/bin/env python
from __future__ import print_function
"""
Different methods to classify genetic mutations.
input: two DNA sequences
output: a list with all mutations of one kind
"""

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
start_table = {"TTG":"M",
               "CTG":"M",
               "ATT":"M",
               "ATC":"M",
               "ATA":"M",
               "ATG":"M",
               "GTG":"M"}

def insertion(x,y):
    mutations = []
    x = x.upper()
    y = y.upper()
    i = 0
    l = 0 # total sequence length (without gaps)
    while i < len(x):
        if x[i] == "-":
            j = i
            has_gap = True
            i += 1
            while has_gap:
                if x[i] != "-":
                    mutations.append("c.{}_{}ins{}".format(l,l+1, y[j:i]))
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
                            "c.{}_{}del{}".format(j+1-k,j+(i-j)-k,x[j:i]))
                    else:
                        mutations.append(
                            "c.{}del{}".format(j+1-k,x[j:i]))
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
                    mutations.append("c.{}{}>{}".format(j+1,x[j:i],y[j:i]))
                    mutated = False
                    i += 1
            if i == len(x):
                mutations.append("c.{}{}>{}".format(j+1,x[j:i],y[j:i]))
        else:
            i += 1
    return(";".join(mutations))

def SNP_coding(x,y,report_silent=False):
    """Return string differences as SNP mutation at protein level"""
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
                mutations.append("p.{}{}_SILENT".format(i+1, a))
            elif b == "*":
                mutations.append("p.{}{}>{}_NONSENSE".format(i+1, a, b))
            elif a != b or a == b == "X":
                mutations.append("p.{}{}>{}_MISSENSE".format(i+1, a, b))
    return(";".join(mutations))

if __name__ == "__main__":
    print("Without silent:",SNP_coding("atgaatggtttgagt","Gtgaatggtttgagt"))
    print("With silent:",SNP_coding("atgaatggtttgagt","Gtgaatggtttgagt",True))
    print("Without silent:",SNP_coding("atgAXGaca","atgAXGact"))
    print("With silent:",SNP_coding("atgAXGaca","atgAXGact",True))
    print(SNP_non_coding("atgAXGaca","atgAXGact"))
    x = "aa---aaaa--aaa-aaa"
    y = "aaACTaaaaGTaaaGaaa"
    print(insertion(x,y))

    x = "aAaa--a--aaaaTGTaa"
    y = "a-aaaaaaaaaaa---aa"
    print(deletion(x,y))
