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

    rf = []
    check = False
    
    for i in range(3):
        splice = seq[i:]

        #removes partial codons by removing bases from end of sequence
        splice = splice[:(len(splice) - len(splice)%3)]
        rf.append(Seq.translate(Seq(splice)))

    for frame in rf:
        if frame[-11] == 'W' and frame[-2:] == 'SS':
            check = True
            break
        else:
            continue

    return check

#sets various nucleotide modification varibles to a number between a specified threshold
#specifically used for deletions
def thresholdSet(l, start, end):

    for i in range(len(l)):
        r = random.randint(start, end)
        if l[i] <= -1:
            l[i] = r
            
    return l

#sets values within threshold if values arent specified within the flag variables
def minMaxDealer(l, start, end):
    check = 0
    for i in range(len(l)):
        if l[i] < 0:
            try:
                r = random.randint(start, end)
            except:
                check = 1
                r = 4
            l[i] = r
    return l, check

#deprecated
'''
#splits a specified insertion gap into 2 seperate values to be handled by each gene seperately
#ie vd insertion = 5 -> v gets 3 nucleotides and d get 2 nucleotides
#randomly assigns a value to both
#specifically used for additions
def insertionDealer(insertion, minI, maxI):

    if insertion < 0:
        r = random.randint(minI, maxI)
        updated = random.randint(0, r)
        first = updated
        second = r - first
        return first, second
    
    else:
        r = random.randint(0, insertion)

    return r, insertion - r
'''

#generate nucleotide sequence based on specificed base parition with <k> bases
def addNuc(k):

    add = [random.randint(0, 99) for i in range(0, k)]
    codonAdd = [basePartition(value, bPoint) for value in add]
    
    return codonAdd

#random nucleotide addition and deletion
#randomizes whether addition or deletion occurs first on the gene
#for addition and deletion to the d gene reverse and unreverse
def randomNucMod(gene, addN, delN):

    rand = random.random()
    addition = ''.join(addNuc(addN))

    #do addition first
    if rand > 0.5:
        modGene = ''.join([gene, addition])[:-delN]
    else:
        modGene = ''.join([gene[:-delN], addition])
        
    return modGene

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
#if additions/deletions are specified, modify the gene respectively
#also specifies the random seed
def rPrint(gList, V, D, J, k, nucVar, m, r):

    #increment variable for sucessful creations of vdj
    inc = 0

    #variable list pertaining to nucleotide modifications
    #ndv = v deletions
    #navd = vd insertions
    #ndd5 = d5 deletions
    #ndd3 = d3 deletions
    #nadj = dj insetions
    #ndj = j deletions
    #m = mutation rate
    
    while inc < k:

        #additions
        addL,checkA = minMaxDealer([nucVar[0], nucVar[1]], nucVar[6], nucVar[7])
        #deletions
        delL,checkD = minMaxDealer([nucVar[2], nucVar[3], nucVar[4], nucVar[5]], nucVar[8], nucVar[9])

        nucVariableList = ', '.join(['ndv:' + str(delL[0]),'navd:' + str(addL[0]),'ndd5:' + str(delL[1]),
                            'ndd3:' + str(delL[2]), 'nadj:' + str(addL[1]),'ndj:' + str(delL[3]), 'm:' + str(m)])

        vdGap = ''.join(addNuc(addL[0]))
        djGap = ''.join(addNuc(addL[1]))

        vKey = random.choice(gList[0])
        dKey = random.choice(gList[1])
        jKey = random.choice(gList[2])

        vGene = ''.join(V[vKey])[:-delL[0]]
        dGene = ''.join(D[dKey])[delL[1]:-delL[2]]
        jGene = ''.join(J[jKey])[delL[3]:]
        
        combStr = "".join([vGene, vdGap, dGene, djGap, jGene])

        #check for stop codons
        combStr = stopChecker(combStr)

        patternCheck = readingFrameChecker(combStr.upper())

        if patternCheck == True:
            inc = inc + 1
        else:
            continue

        if checkA == 1:
            print ('* invalid start/end conditions for addition interval, defaulting to 4')
        if checkD == 1:
            print ('* invalid start/end conditions for deletion interval, defaulting to 4')
            
        print ('>' + ", ".join(('v:' + vKey, 'd:' + dKey, 'j:' + jKey, nucVariableList)))
        mutStr = mutation(combStr, m)

        print (mutStr)

    print ('\nrandom seed: ' + str(r))

    
def main():
    parser = argparse.ArgumentParser(description="Generate VDJ recombinations with junctional diversity "
                                     "and somatic hypermutations. Tests are in place to ensure that "
                                     "sequences do not contain stop codons and that J gene is in the "
                                     "correct reading frame. Default values for optional command "
                                     "line argument are list in parentheses below. "
                                     "Each sequence is annotated with: V allele, D allele, J allele, "
                                     "V deletions, VD insertions, D5' deletions, D3' deletions, "
                                     "DJ insertions, J deletions")

    #vdj files args
    parser.add_argument('-v', '--vgene', dest='vGene', type = validFile, required = True,
                        help = "V gene .fasta file")
    parser.add_argument('-d', '--dgene', dest='dGene', type = validFile, required = True,
                        help = "D gene .fasta file")
    parser.add_argument('-j', '--jgene', dest='jGene', type = validFile, required = True,
                        help = "J gene .fasta file")
    parser.add_argument('-a', '--allele', dest='aFile', required = False, help = "Total allele .txt file"
                        " and their respective allele. Currently formatted as just the gene with the "
                        " attached allele followed by new line character. Ex: IGHD1*01\\nIGHD2*01\\n...")

    #nucleotide modification args
    parser.add_argument('-maxi', dest = 'maxI', type = int, default = 4,
                        help = "max # insertions possible for each junction (4)")
    parser.add_argument('-mini', dest = 'minI', type = int, default = 4,
                        help = "min # insertions possible for each junction (4)")
    parser.add_argument('-maxd', dest = 'maxD', type = int, default = 4,
                        help = "max # deletions possible for each junction (4)")
    parser.add_argument('-mind', dest = 'minD', type = int, default = 4,
                        help = "min # deletions possible for each junction (4)")


    parser.add_argument('-ndv', dest='ndv', type = int, default = -1,
                        help = "# nt to be deleted from the V gene. Set >= 0 to override maxd/mind (-1)")
    parser.add_argument('-navd', dest='navd', type = int, default = -1,
                        help = "# nt to be added to the VD junction. Set >= 0 to override maxi/mini (-1)")
    parser.add_argument('-ndd5', dest='ndd5', type = int, default = -1,
                        help = "# nt to be deleted from 5' end D gene. Set >= 0 to override maxd/mind (-1)")
    parser.add_argument('-ndd3', dest='ndd3', type = int, default = -1,
                        help = "# nt to be deleted the 3' end D gene. Set >= 0 to override maxd/mind (-1)")
    parser.add_argument('-nadj', dest='nadj', type = int, default = -1,
                        help = "# nt to be added to the DJ junction. Set >= 0 to override maxi/mini (-1)")
    parser.add_argument('-ndj', dest='ndj', type = int, default = -1,
                        help = "# nt to be deleted from the J gene. Set >= 0 to override maxd/mind (-1)")
    
    #misc arg variables
    parser.add_argument('-m', dest = 'mut', type = int, default = 0,
                        help = "pct of each base mutating into another base (0)")
    parser.add_argument('-r', dest='rSeed', type = int, help = "Optional random seed input"
                        " Otherwise will generate seed for the current data")
    parser.add_argument('-k', dest='kIter', type = int, required = True, help = "# of iterations")
    
    args = parser.parse_args()

    #Checks for random seed
    if args.rSeed:
        rNum = args.rSeed
    else:
        dt = datetime.now()
        rNum = dt.microsecond

    #applys random seed to every possible random var
    random.seed(rNum)
        
    m = percentMut(args.mut)
    
    vFile = open(args.vGene, 'r')
    dFile = open(args.dGene, 'r')
    jFile = open(args.jGene, 'r')

    dictV = parse(vFile)
    dictD = parse(dFile)
    dictJ = parse(jFile)

    if args.aFile:
        aFile = open(args.aFile, 'r')
        gList = compare(aFile)
    else:
        gList = (list(dictV.keys()), list(dictD.keys()), list(dictJ.keys()))

    nucVariables = [args.navd, args.nadj, args.ndv, args.ndd5,
            args.ndd3, args.ndj, args.minI, args.maxI, args.minD, args.maxD]

    rPrint (gList, dictV, dictD, dictJ, args.kIter, nucVariables, m, rNum)


main()
