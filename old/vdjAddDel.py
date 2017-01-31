import sys
import os
import argparse
import random
from datetime import datetime
from bisect import bisect


#checks for valid files within argparse V, D, J selection
def validFile (fileName):
    base, ext = os.path.splitext(fileName)
    if ext.lower() not in ('.fasta'):
        raise argparse.ArgumentTypeError('Error, requires .fasta file')

    return fileName

#prevents frame shift mutations by ensuring that the user inputted values
#for both additions and deletions have a (delta) that is divisible by 3
def frameBlocker(a1, d1, a2, d2):
    
    if ((a1-d1)%3 + (a2-d2)%3 == 0):
        return 1

    while ((a1-d1)%3 != 0):
        print ('Enter new values for Add1 and Del1 so that their difference '
               'when modded by 3 produces a result of 0')
        a1 = int(input ('Enter a new value for add1: '))
        d1 = int(input ('Enter a new value for del1: '))

    while ((a2-d2)%3 != 0):
        print ('Enter new values for Add2 and Del2 so that their difference '
               'when modded by 3 produces a result of 0')
        a2 = int(input ('Enter a new value for add2: '))
        d2 = int(input ('Enter a new value for del2: '))
    
    return 1       
            
    

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
#if additions/deletions are specified, compute and add those in front or behind
#the d gene when printing
#also specifies the random seed
def rPrint(gList, V, D, J, k, a1, d1, a2, d2, r):

    bPoint = [25, 50, 75, 100]

    if isinstance(r, bool):
        dt = datetime.now()
        r = dt.microsecond

    random.seed(r)
    
    for i in range(0, k):
            
        vdBP = [random.randint(0, 99) for i in range(0, a1)]
        djBP = [random.randint(0, 99) for i in range(0, a2)]
    
        vdAdd = [basePartition(value, bPoint) for value in vdBP]
        djAdd = [basePartition(value, bPoint) for value in djBP]
        
        vKey = random.choice(gList[0])
        dKey = random.choice(gList[1])
        jKey = random.choice(gList[2])

    
        vdAdd = ''.join(vdAdd)
        djAdd = ''.join (djAdd)
        
        print ('>' + ", ".join((vKey, dKey, jKey)))
        unchangedD = D[dKey]
        delD = unchangedD[d1:-d2]
        print ("".join((V[vKey], vdAdd, delD, djAdd, J[jKey])))

    print ('\nrandom seed: ' + str(r))

    
def main():
    parser = argparse.ArgumentParser(description='Process specific alleles of V, D, and J')
    
    parser.add_argument('-v', '--vgene', dest='vGene', type = validFile, required = True, help = "V gene .fasta file")
    parser.add_argument('-d', '--dgene', dest='dGene', type = validFile, required = True, help = "D gene .fasta file")
    parser.add_argument('-j', '--jgene', dest='jGene', type = validFile, required = True, help = "J gene .fasta file")

    parser.add_argument('-k', dest='kIter', type = int, required = True, help = "# of iterations")

    parser.add_argument('-a', '--allele', dest='aFile', required = True, help = "Total allele .txt file"
                        " and their respective allele. Currently formatted as just the gene with the "
                        " attached allele followed by new line character. Ex: IGHD1*01\\nIGHD2*01\\n...")

    parser.add_argument('-na1', dest='nAdd1', type = int, default = 0,
                        help = "# of nucleotides to be added inbetween V and D gene, default = 0")
    parser.add_argument('-nd1', dest='nDel1', type = int, default = 0,
                        help = "# of nucleotides to be deleted inbetween V and D gene, default = 0")
    parser.add_argument('-na2', dest='nAdd2', type = int, default = 0,
                        help = "# of nucleotides to be added inbetween D and J gene, default = 0")
    parser.add_argument('-nd2', dest='nDel2', type = int, default = 0,
                        help = "# of nucleotides to be deleted inbetween D and J gene, default = 0")

    parser.add_argument('-r', dest='rSeed', type = int, help = "Optional random seed input"
                        " Otherwise will generate seed for the current data")
    
    args = parser.parse_args()

    #Checks for random seed
    if args.rSeed:
        rNum = args.rSeed
    else:
        rNum = False
        
    frameBlocker(args.nAdd1, args.nDel1, args.nAdd2, args.nDel2)

    aFile = open(args.aFile, 'r')
    vFile = open(args.vGene, 'r')
    dFile = open(args.dGene, 'r')
    jFile = open(args.jGene, 'r')

    dictV = parse(vFile)
    dictD = parse(dFile)
    dictJ = parse(jFile)
    
    gList = compare(aFile)

    rPrint (gList, dictV, dictD, dictJ, args.kIter, args.nAdd1, args.nDel1,
            args.nAdd2, args.nDel2, rNum)


main()
