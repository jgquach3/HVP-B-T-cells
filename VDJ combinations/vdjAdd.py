import sys
import os
import argparse
import random
import math
from Bio.Seq import Seq
from datetime import datetime
from bisect import bisect


#Globals#

#percentages of each base pair, can be changed to simulate different percentages
#current order = ACTG
bPoint = [25, 50, 75, 100]

#checks for valid files within argparse V, D, J selection
def validFile (fileName):
    base, ext = os.path.splitext(fileName)
    if ext.lower() not in ('.fasta'):
        raise argparse.ArgumentTypeError('Error, requires .fasta file')

    return fileName
  

#checks to see if the user input of the mutation rate is valid: between 0 - 100%
def percentMut (m):

    while (m > 100 or m < 0):
        print ('Mutation rate is a percentage between 0 and 100')

        m = int(input('Enter a new value for the mutation rate: '))

    return int(m)

#checks for stop codons and randomly changes them to another protein
def stopChecker(chain):

    #arbitrarily large int to handle multiple creations of stop codons
    i = 1000000000
    
    check = True
    while check == True:

        if 'tag' in chain:

            nC = [random.randint(0, 99) for i in range(0, 3)]
            nCodon = [basePartition(value, bPoint) for value in nC]
            nCodon = ''.join(nCodon)
            chain = chain.replace('tag', nCodon, i)
            continue
                
        elif 'taa' in chain:

            nC = [random.randint(0, 99) for i in range(0, 3)]
            nCodon = [basePartition(value, bPoint) for value in nC]
            nCodon = ''.join(nCodon)
            chain = chain.replace('taa', nCodon, i)
            continue
        
        elif 'tga' in chain:

            nC = [random.randint(0, 99) for i in range(0, 3)]
            nCodon = [basePartition(value, bPoint) for value in nC]
            nCodon = ''.join(nCodon)
            chain = chain.replace('tga', nCodon, i)
            continue

        else:
            
            check = False

    return chain
    

#partitions a base pair list based on set percentages
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


#checks to see if sequence belongs to reading frame:
#1 - no additional nucleotides
#2 - 1 additional nucleotide
#3 - 2 additional nucleotides
#and if it satisfies the following pattern: {X} - {W} - {X}8 - {SS}
def readingFrameChecker(seq):

    rf1 = Seq.translate(Seq(seq))
    rf2 = Seq.translate(Seq(str(seq[1:])))
    rf3 = Seq.translate(Seq(str(seq[2:])))

    check = False

    for frame in [rf1, rf2, rf3]:

        if frame[-11] == 'W' and frame[-2:] == 'SS':

            check = True
            break
        else:
            continue

    return check

#generate nucleotide sequence based on specificed base parition with <k> bases
def addNuc(k):

    add = [random.randint(0, 99) for i in range(0, k)]
    codonAdd = [basePartition(value, bPoint) for value in add]
    
    return codonAdd


#applies the mutation percentage to every base within the string and
#mutates the bases if applicable
#also checks if there are any stop codons and adjusts them if there are
def mutation (chain, m):

    chain = list(chain)
    
    for i in range(len(chain)):

        rate = random.randint(0, 101)

        if (rate <= m):

            baseNumber = random.randint(0,99)
            nBase = basePartition(baseNumber, bPoint)
            chain[i] = nBase

    #check for stop codons
    chain = stopChecker(''.join(chain))

    return chain

        
#randomly takes 1 value from each V,D, and J list and generates sequences
#based on k iterations
#if additions/deletions are specified, compute and add those in front or behind
#the d gene when printing
#also specifies the random seed
def rPrint(gList, V, D, J, k, a1, d1, a2, d2, m, r):

    if isinstance(r, bool):
        dt = datetime.now()
        r = dt.microsecond

    random.seed(r)

    #increment variable for sucessful creations of vdj
    inc = 0
    
    while inc < k:
    
        vdAdd = addNuc(a1)
        djAdd = addNuc(a2)
        
        vKey = random.choice(gList[0])
        dKey = random.choice(gList[1])
        jKey = random.choice(gList[2])

    
        vdAdd = ''.join(vdAdd)
        djAdd = ''.join(djAdd)
        
        print ('>' + ", ".join((vKey, dKey, jKey)))

        #only deleted on D gene, does not address possible deletion on V/J
        unchangedD = D[dKey]
        if d1 == 0 and d2 == 0:
            delD = unchangedD
        else:
            delD = unchangedD[d1:-d2]

        combStr = "".join((V[vKey], vdAdd, delD, djAdd, J[jKey]))

        #check for stop codons
        combStr = stopChecker(combStr)

        patternCheck = readingFrameChecker(combStr.upper())

        if patternCheck == True:
            inc = inc + 1
        else:
            continue
        
        mutStr = mutation(combStr, m)        
        print (mutStr)

    print ('\nrandom seed: ' + str(r))

    
def main():
    parser = argparse.ArgumentParser(description='Process specific alleles of V, D, and J')
    
    parser.add_argument('-v', '--vgene', dest='vGene', type = validFile, required = True,
                        help = "V gene .fasta file")
    parser.add_argument('-d', '--dgene', dest='dGene', type = validFile, required = True,
                        help = "D gene .fasta file")
    parser.add_argument('-j', '--jgene', dest='jGene', type = validFile, required = True,
                        help = "J gene .fasta file")

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

    parser.add_argument('-m', dest = 'mut', type = int, default = 0,
                        help = "% of each base mutating into another base")

    parser.add_argument('-r', dest='rSeed', type = int, help = "Optional random seed input"
                        " Otherwise will generate seed for the current data")
    
    args = parser.parse_args()

    #Checks for random seed
    if args.rSeed:
        rNum = args.rSeed
    else:
        rNum = False
        
    m = percentMut(args.mut)

    aFile = open(args.aFile, 'r')
    vFile = open(args.vGene, 'r')
    dFile = open(args.dGene, 'r')
    jFile = open(args.jGene, 'r')

    dictV = parse(vFile)
    dictD = parse(dFile)
    dictJ = parse(jFile)
    
    gList = compare(aFile)

    rPrint (gList, dictV, dictD, dictJ, args.kIter, args.nAdd1, args.nDel1,
            args.nAdd2, args.nDel2, m, rNum)


main()
