import sys
import os
import argparse
import csv
import random
import math

cdr = {}

#construct unique CDR dictionary using the CDR from each specified file
#increment counts CDR is encountered again
def populate(f, cdr):

    for line in f:
        element = line[1]
        
        if element in cdr:
            value = cdr[element]
            cdr[element] = value + 1
        elif len(element) < 1:
            pass
        else:
            cdr[element] = 1

    return

def main():
    parser = argparse.ArgumentParser(description='Count number of CDRs (Adaptive Seq)')

    parser.add_argument('-f', '--file', action='append', dest='files', default=[],
                        help='Add files to list')

    args = parser.parse_args()

    files = args.files

    for i in range(len(files)):
        file = files[i]

        f = csv.reader(open(file, 'rb'), delimiter='\t')
        
        populate(f, cdr)

    output = sorted(cdr.items(), reverse = True, key=lambda x:x[1])
    o = open('out.txt', 'w')
    
    for item in output:
        o.write(item[0] + ' ' + str(item[1]) + '\n')
    
    
main()
