import sys
import os
import argparse
import csv
import random
import math

#global variables for the dictionaries
codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
cdr = {}

#construct unique CDR dictionary using the CDR from each specified file
def populate(f, cdr):

    for item in f[1:]:
        element = item[1]
        
        if element in cdr:
            pass
        elif len(element) < 1:
            pass
        else:
            cdr[element] = 0

    return

#convert seq to its protein form for comparison with CDR
def translate(text):

    l = []
    
    for i in range(0, len(text), 3):
        codon = text[i:i+3]
        letter = codontable[codon]
        l.append(letter)

    return ''.join(l)

#utilizing the populated CDR dictionary, incrememnt the count of
#a CDR if its present within a sequence
def compare(f, cdr):

    k = cdr.keys()
    
    for item in f[1:]:
        seq = translate(item[0])
        
        for key in k:
            if key in seq:
                value = cdr[key]
                cdr[key] = value + 1

    return
        

def main():
    parser = argparse.ArgumentParser(description='Count number of CDRs (Adaptive Seq)')

    parser.add_argument('-f', '--file', action='append', dest='files', default=[],
                        help='Add files to list')

    args = parser.parse_args()

    files = args.files

    for i in range(len(files)):
        file = files[i]

        f = list(csv.reader(open(file, 'rb'), delimiter='\t'))
        
        populate(f, cdr)
        compare(f, cdr)

    output = sorted(cdr.items(), reverse = True, key=lambda x:x[1])
    o = open('out.txt', 'w')
    
    for item in output:
        o.write(item[0] + ' ' + str(item[1]) + '\n')
    
    
main()
