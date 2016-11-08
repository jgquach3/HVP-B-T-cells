import sys
import os
import argparse
from datetime import datetime
import random

#checks for valid files within argparse V, D, J selection
def validFile (fileName):
    base, ext = os.path.splitext(fileName)
    if ext.lower() not in ('.fasta'):
        raise argparse.ArgumentTypeError('Error, requires .fasta file')

    return fileName

#parses V, D, J files and concats their base pairs together
def parse(fileP):

    dictP = {}
    key = ''
    
    for line in fileP:
        
        if line.startswith('>'):
            key = line.lstrip('>').rstrip('\n')
            dictP.update({key: ''})
        else:
            line = line.rstrip('\n')
            dictP.update({key: "".join((dictP[key], line))})

    return dictP

#creates separate lists of V/D/J genes depending on allele present
def compare(aFile):

    combine = ''
    vList = []
    dList = []
    jList = []
    
    for line in aFile:
        
        selector = line[3]
        line = line.rsplit('\n')

        if selector == 'V':
            vList.append(line[0])
        elif selector == 'D':
            dList.append(line[0])
        else:
            jList.append(line[0])

    return vList, dList, jList

#randomly takes 1 value from each V,D, and J list and generates sequences
#based on k iterations, also creates a seed based on the current time
#in microseconds if no seed is provided

#output: > V key, D key, J key, | random seed
def rPrint(gList, V, D, J, k, r):
    
    if isinstance(r, bool):
        dt = datetime.now()
        r = dt.microsecond

    random.seed(r)

    for i in range(0, k):
        vKey = random.choice(gList[0])
        dKey = random.choice(gList[1])
        jKey = random.choice(gList[2])

        print ('>' + ", ".join((vKey, dKey, jKey)) + " | " + str(r))
        print ("".join((V[vKey], D[dKey], J[jKey])))

#takes 6 arguments
#3 .fasta files for the vgenes, dgenes, and j genes
#1 .txt file that contains a subset of genes and their alleles
#1 int that determines how many iterations of sequences
#1 int the random seed
        
def main():
    parser = argparse.ArgumentParser(description='Process specific alleles of V, D, and J')
    
    parser.add_argument('-v', '--vgene', dest='vGene', type = validFile, help = "V gene .fasta file")
    parser.add_argument('-d', '--dgene', dest='dGene', type = validFile, help = "D gene .fasta file")
    parser.add_argument('-j', '--jgene', dest='jGene', type = validFile, help = "J gene .fasta file")

    parser.add_argument('-a', '--allele', dest='aFile', help = "Total allele .txt file"
                        " and their respective allele. Currently formatted as just the gene with the "
                        " attached allele followed by new line character. Ex: IGHD1*01\\nIGHD2*01\\n...")

    parser.add_argument('-k', dest='kIter', type = int, help = "# of iterations")
    
    parser.add_argument('-r', dest='rSeed', type = int, help = "Optional random seed input"
                        " Otherwise will generate seed for the current data")


    args = parser.parse_args()

    if args.rSeed:
        rNum = args.rSeed
    else:
        rNum = False

    vFile = open(args.vGene, 'r')
    dFile = open(args.dGene, 'r')
    jFile = open(args.jGene, 'r')
    aFile = open(args.aFile, 'r')

    dictV = parse(vFile)
    dictD = parse(dFile)
    dictJ = parse(jFile)

    gList = compare(aFile)

    rPrint (gList, dictV, dictD, dictJ, args.kIter, rNum)


main()
