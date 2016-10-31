import sys
import os
import argparse
import random


#checks for valid files within argparse
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
            

#randomly select 1 key from each dict, outputting the key and concatted values            
def randomCombine(V, D, J, k):      

    for i in range(0, k):
        vKey = random.choice(list(V.keys()))
        dKey = random.choice(list(D.keys()))
        jKey = random.choice(list(J.keys()))

        print ('>' + ", ".join((vKey, dKey, jKey)))
        print ("".join((V[vKey], D[dKey], J[jKey])))
    
def main():

    #defining the argument parser
    parser = argparse.ArgumentParser(description='Process V, D, and J files with k iterations')
    
    parser.add_argument('-v', '--vgene', dest='vGene', type = validFile, help = "V gene .fasta file here")
    parser.add_argument('-d', '--dgene', dest='dGene', type = validFile, help = "D gene .fasta file here")
    parser.add_argument('-j', '--jgene', dest='jGene', type = validFile, help = "J gene .fasta file here")
    parser.add_argument('-k', dest='kIter', type = float, help = "# of iterations")
    
                                 
    args = parser.parse_args()


    fileV = open(args.vGene, 'r')
    fileD = open(args.dGene, 'r')
    fileJ = open(args.jGene, 'r')

    dictV = parse(fileV)
    dictD = parse(fileD)
    dictJ = parse(fileJ)

    randomCombine(dictV, dictD, dictJ, int(args.kIter))

main()
