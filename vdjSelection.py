import sys
import os
import argparse
import random
from bisect import bisect


#checks for valid files within argparse V, D, J selection
def validFile (fileName):
    base, ext = os.path.splitext(fileName)
    if ext.lower() not in ('.fasta'):
        raise argparse.ArgumentTypeError('Error, requires .fasta file')

    return fileName

#partitions a bp list based on set percentages
def basePartition (value, bPoint):

    bp = 'atcg'
    i = bisect(bPoint, value)

    return bp[i]

    
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
#based on k iterations
#if additions are specified, compute and add those in between as well
def rPrint(gList, V, D, J, k, t1, t2):

    bPoint = [25, 50, 75, 100]
    
    for i in range(0, k):

        vdBP = [random.randint(0, 99) for i in range(0, t1*3)]
        djBP = [random.randint(0, 99) for i in range(0, t2*3)]
    
        vdAdd = ''.join([basePartition(value, bPoint) for value in vdBP])
        djAdd = ''.join([basePartition(value, bPoint) for value in djBP])
        
        vKey = random.choice(gList[0])
        dKey = random.choice(gList[1])
        jKey = random.choice(gList[2])

        print ('>' + ", ".join((vKey, dKey, jKey)))
        print ("".join((V[vKey], vdAdd, D[dKey], djAdd, J[jKey])))

    
def main():
    parser = argparse.ArgumentParser(description='Process specific alleles of V, D, and J')
    
    parser.add_argument('-v', '--vgene', dest='vGene', type = validFile, help = "V gene .fasta file")
    parser.add_argument('-d', '--dgene', dest='dGene', type = validFile, help = "D gene .fasta file")
    parser.add_argument('-j', '--jgene', dest='jGene', type = validFile, help = "J gene .fasta file")

    parser.add_argument('-k', dest='kIter', type = int, help = "# of iterations")

    parser.add_argument('-a', '--allele', dest='aFile', help = "Total allele .txt file"
                        " and their respective allele. Currently formatted as just the gene with the "
                        " attached allele followed by new line character. Ex: IGHD1*01\\nIGHD2*01\\n...")

    parser.add_argument('-tAdd1', dest='tAdd1', type = int, default = 0,
                        help = "# of triplets to be added inbetween V and D gene, default = 0")
    parser.add_argument('-tAdd2', dest='tAdd2', type = int, default = 0,
                        help = "# of triplets to be added inbetween D and J gene, default = 0")
    
    args = parser.parse_args()


    aFile = open(args.aFile, 'r')
    vFile = open(args.vGene, 'r')
    dFile = open(args.dGene, 'r')
    jFile = open(args.jGene, 'r')

    dictV = parse(vFile)
    dictD = parse(dFile)
    dictJ = parse(jFile)
    
    gList = compare(aFile)

    rPrint (gList, dictV, dictD, dictJ, args.kIter, args.tAdd1, args.tAdd2)


main()
